use monstertruck_geometry::prelude::*;
use rustc_hash::{FxHashMap as HashMap, FxHashSet as HashSet};
use smallvec::SmallVec;

use super::error::DraftError;
use super::types::*;

type Result<T> = std::result::Result<T, DraftError>;

/// A plane equation `normal . x = distance`, in absolute (world)
/// coordinates. `normal` is always a unit vector.
type PlaneEq = (Vector3, f64);

/// A face's replacement plane: the plane equation used to reconcile shared
/// vertices, plus a point known to lie on that plane, used to rebuild the
/// face's own surface geometry.
#[derive(Clone, Copy)]
struct FacePlane {
    normal: Vector3,
    distance: f64,
    point: Point3,
}

/// Finds the point satisfying every plane equation in `constraints`, using
/// `reference` to resolve the under-determined cases (1 or 2 constraints).
/// This is the same plane-intersection problem the shell/thicken module
/// solves for its own vertex reconciliation, generalized from "offset every
/// plane by a scalar along its own normal" to "replace each plane with an
/// arbitrary new one": here `constraints` carries each adjacent face's
/// possibly-rotated plane directly, in absolute `(normal, distance)` form,
/// rather than a per-face scalar offset.
///
/// - **1 constraint**: the orthogonal projection of `reference` onto the
///   plane.
/// - **2 constraints**: the point reached from `reference` by a displacement
///   confined to `span(normal0, normal1)` -- perpendicular to the planes'
///   shared line, so a `reference` that lies on both planes' original shared
///   edge keeps that edge's direction.
/// - **3 constraints**: the unique intersection point of the three planes
///   (independent of `reference`).
fn intersect_planes(reference: Point3, constraints: &[PlaneEq]) -> Result<Point3> {
    match *constraints {
        [(normal, distance)] => {
            let excess = normal.dot(reference.to_vec()) - distance;
            Ok(reference - excess * normal)
        }
        [(normal0, distance0), (normal1, distance1)] => {
            let cos = normal0.dot(normal1);
            let det = 1.0 - cos * cos;
            if det.abs() < TOLERANCE2 {
                Err(DraftError::DegenerateVertex)
            } else {
                let local0 = distance0 - normal0.dot(reference.to_vec());
                let local1 = distance1 - normal1.dot(reference.to_vec());
                let a = (local0 - cos * local1) / det;
                let b = (local1 - cos * local0) / det;
                Ok(reference + a * normal0 + b * normal1)
            }
        }
        [
            (normal0, distance0),
            (normal1, distance1),
            (normal2, distance2),
        ] => {
            let matrix = Matrix3::from_cols(normal0, normal1, normal2).transpose();
            if matrix.determinant().abs() < TOLERANCE {
                Err(DraftError::DegenerateVertex)
            } else {
                // SAFETY: the determinant was just checked to be non-zero.
                let inverse = matrix.invert().unwrap();
                Ok(Point3::from_vec(
                    inverse * Vector3::new(distance0, distance1, distance2),
                ))
            }
        }
        _ => Err(DraftError::UnsupportedVertexDegree(constraints.len())),
    }
}

/// Collects the distinct plane constraints of a vertex, merging near-equal
/// normals and rejecting coplanar faces that request conflicting planes.
fn dedup_constraints(raw: impl IntoIterator<Item = PlaneEq>) -> Result<SmallVec<[PlaneEq; 4]>> {
    let mut constraints: SmallVec<[PlaneEq; 4]> = SmallVec::new();
    for (normal, distance) in raw {
        match constraints
            .iter()
            .find(|&&(other, _)| (other - normal).magnitude2() < TOLERANCE2)
        {
            Some(&(_, other_distance)) if (other_distance - distance).abs() < TOLERANCE => {}
            Some(_) => return Err(DraftError::DegenerateVertex),
            None => constraints.push((normal, distance)),
        }
    }
    Ok(constraints)
}

/// Builds a plane through `point` with the given unit `normal`, choosing an
/// arbitrary but consistent in-plane basis.
fn plane_through(normal: Vector3, point: Point3) -> Plane {
    let helper = match normal.x.abs() < 0.9 {
        true => Vector3::unit_x(),
        false => Vector3::unit_y(),
    };
    let u = normal.cross(helper).normalize();
    // `u` is perpendicular to `normal`, so `u x v = normal` exactly.
    let v = normal.cross(u);
    Plane::new(point, point + u, point + v)
}

/// The oriented replacement plane for one face: unchanged if `drafted` is
/// `false`, otherwise the face's original plane rotated by `angle` about its
/// intersection line with the neutral plane, tilting so it opens toward
/// `pull_direction`.
fn face_plane(
    face: &Face,
    drafted: bool,
    neutral_normal: Vector3,
    neutral_distance: f64,
    pull_direction: Vector3,
    angle: Rad<f64>,
) -> Result<FacePlane> {
    let normal = face.oriented_surface().normal();
    let origin = face.surface().origin();
    let distance = normal.dot(origin.to_vec());
    if !drafted {
        return Ok(FacePlane {
            normal,
            distance,
            point: origin,
        });
    }
    let axis_raw = normal.cross(neutral_normal);
    if axis_raw.magnitude() < TOLERANCE {
        return Err(DraftError::ParallelToNeutralPlane);
    }
    // The hinge line is `face plane ∩ neutral plane`; this point on it is
    // fixed by the rotation and lets us both check for the degenerate case
    // above through the shared solver and place the rotated plane exactly.
    let hinge_point = intersect_planes(
        Point3::origin(),
        &[(normal, distance), (neutral_normal, neutral_distance)],
    )
    .map_err(|_| DraftError::ParallelToNeutralPlane)?;
    let axis = axis_raw.normalize();
    // `pull_direction`'s component along `axis x normal` says which way a
    // point far from the hinge line, in the pull direction, moves as the
    // plane rotates; flip the rotation to make that motion inward.
    let alignment = pull_direction.dot(axis.cross(normal));
    if alignment.abs() < TOLERANCE {
        return Err(DraftError::AmbiguousPullDirection);
    }
    let signed_angle = Rad(angle.0 * alignment.signum());
    let rotation = Matrix3::from_axis_angle(axis, signed_angle);
    let new_normal = (rotation * normal).normalize();
    let new_distance = new_normal.dot(hinge_point.to_vec());
    Ok(FacePlane {
        normal: new_normal,
        distance: new_distance,
        point: hinge_point,
    })
}

/// Rebuilds `shell` with the faces in `drafted_faces` tilted by `angle`
/// about their intersection with `neutral_plane`, opening toward
/// `pull_direction`. Every vertex and edge is reconciled against the
/// (possibly unchanged) planes of its adjacent faces.
pub(super) fn drafted_shell(
    shell: &Shell,
    drafted_faces: &HashSet<usize>,
    neutral_plane: Plane,
    pull_direction: Vector3,
    angle: Rad<f64>,
) -> Result<Shell> {
    let neutral_normal = neutral_plane.normal();
    let neutral_distance = neutral_normal.dot(neutral_plane.origin().to_vec());

    let face_planes = shell
        .iter()
        .enumerate()
        .map(|(face_index, face)| {
            face_plane(
                face,
                drafted_faces.contains(&face_index),
                neutral_normal,
                neutral_distance,
                pull_direction,
                angle,
            )
        })
        .collect::<Result<Vec<_>>>()?;

    let mut vertex_faces: HashMap<VertexId, (Vertex, SmallVec<[usize; 4]>)> = HashMap::default();
    shell.iter().enumerate().for_each(|(face_index, face)| {
        face.vertex_iter().for_each(|vertex| {
            let entry = vertex_faces
                .entry(vertex.id())
                .or_insert_with(|| (vertex.clone(), SmallVec::new()));
            if !entry.1.contains(&face_index) {
                entry.1.push(face_index);
            }
        });
    });

    let vertices = vertex_faces
        .into_iter()
        .map(|(vertex_id, (vertex, face_indices))| {
            let raw = face_indices
                .iter()
                .map(|&index| (face_planes[index].normal, face_planes[index].distance));
            let constraints = dedup_constraints(raw)?;
            let new_point = intersect_planes(vertex.point(), &constraints)?;
            Ok((vertex_id, Vertex::new(new_point)))
        })
        .collect::<Result<HashMap<_, _>>>()?;

    let mut edges: HashMap<EdgeId, Edge> = HashMap::default();
    for edge in shell.edge_iter() {
        if edges.contains_key(&edge.id()) {
            continue;
        }
        let (front, back) = edge.absolute_ends();
        let lookup = |vertex: &Vertex| {
            vertices
                .get(&vertex.id())
                .ok_or(DraftError::GeometryFailed {
                    context: "drafted vertex lookup for edge",
                })
        };
        let (new_front, new_back) = (lookup(front)?, lookup(back)?);
        let original = back.point() - front.point();
        let drafted_direction = new_back.point() - new_front.point();
        if drafted_direction.magnitude() < TOLERANCE || drafted_direction.dot(original) <= 0.0 {
            return Err(DraftError::ExcessiveAngle { radians: angle.0 });
        }
        let line = Line(new_front.point(), new_back.point());
        edges.insert(edge.id(), Edge::new(new_front, new_back, line));
    }

    let faces = shell
        .iter()
        .enumerate()
        .map(|(face_index, face)| {
            // `face.boundaries()` (not `absolute_boundaries()`) is the wire
            // already corrected for this face's orientation, so it is always
            // outward-consistent with `face.oriented_surface()` -- the same
            // convention `face_plane` used to compute `plane_info.normal`
            // (via `face.oriented_surface().normal()`, rotated). `plane`
            // below is built by `plane_through` from that already-outward
            // normal, so pairing it with this wire and leaving the new
            // face's orientation at its `Face::try_new` default of `true`
            // keeps both members outward-consistent regardless of whatever
            // orientation the *original* face happened to have; re-applying
            // `face.orientation()` here (as a previous version of this code
            // did, pairing the new plane with the *unflipped*
            // `absolute_boundaries()` instead) would flip only the wire and
            // not the (already correctly outward) surface for exactly the
            // faces whose original orientation was `false`, silently
            // inverting the rebuilt face's outward direction.
            let boundaries = face
                .boundaries()
                .iter()
                .map(|wire| {
                    wire.iter()
                        .map(|edge| {
                            let new_edge =
                                edges.get(&edge.id()).ok_or(DraftError::GeometryFailed {
                                    context: "drafted edge lookup for face",
                                })?;
                            Ok(match edge.orientation() {
                                true => new_edge.clone(),
                                false => new_edge.inverse(),
                            })
                        })
                        .collect::<Result<Wire>>()
                })
                .collect::<Result<Vec<_>>>()?;
            let plane_info = &face_planes[face_index];
            let plane = plane_through(plane_info.normal, plane_info.point);
            let new_face = Face::try_new(boundaries, plane)
                .map_err(|source| DraftError::InvalidOutputTopology { source })?;
            Ok(new_face)
        })
        .collect::<Result<Vec<_>>>()?;

    Ok(faces.into_iter().collect())
}
