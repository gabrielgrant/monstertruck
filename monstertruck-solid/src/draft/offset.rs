use monstertruck_geometry::prelude::*;
use rustc_hash::{FxHashMap as HashMap, FxHashSet as HashSet};
use smallvec::SmallVec;

use crate::planar::{dedup_constraints, intersect_planes};

use super::error::DraftError;
use super::types::*;

type Result<T> = std::result::Result<T, DraftError>;

/// A face's replacement plane: the plane equation used to reconcile shared
/// vertices, plus a point known to lie on that plane, used to rebuild the
/// face's own surface geometry.
#[derive(Clone, Copy)]
struct FacePlane {
    normal: Vector3,
    distance: f64,
    point: Point3,
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
