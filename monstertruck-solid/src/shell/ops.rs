use monstertruck_geometry::prelude::*;
use rustc_hash::{FxHashMap as HashMap, FxHashSet as HashSet};

use super::error::ShellError;
use super::offset::{OffsetElements, offset_elements};
use super::types::*;

type Result<T> = std::result::Result<T, ShellError>;
type ExternalShell<C, S> = monstertruck_topology::Shell<Point3, C, S>;
type ExternalSolid<C, S> = monstertruck_topology::Solid<Point3, C, S>;
type ExternalFaceId<S> = monstertruck_topology::FaceId<S>;

/// Thickens a planar-faced shell by a uniform `thickness` into a solid.
///
/// Every face is offset by `thickness` against its oriented normal, so the
/// original shell remains one boundary of the material and the offset shell
/// becomes the other. For an open shell, planar side walls are built along
/// the boundary wires; for a closed shell, the result is a hollow solid
/// with two boundary shells (equivalent to [`hollow`] with no open faces).
///
/// # Errors
/// - [`ShellError::UnsupportedGeometry`]: a surface is not a
///   [`Plane`](monstertruck_geometry::prelude::Plane), an edge curve is not
///   a [`Line`](monstertruck_geometry::prelude::Line), or a boundary vertex
///   offsets out of its side-wall plane.
/// - [`ShellError::NotOrientable`]: the shell is disconnected, non-manifold,
///   or not orientable.
/// - [`ShellError::NonPositiveThickness`] / [`ShellError::ThicknessTooLarge`]:
///   invalid thickness.
/// - [`ShellError::DegenerateVertex`] / [`ShellError::UnsupportedVertexDegree`]:
///   a vertex offset cannot be reconciled by plane intersection.
pub fn thicken<C: ThickenableCurve, S: ThickenableSurface>(
    shell: &ExternalShell<C, S>,
    thickness: f64,
) -> Result<ExternalSolid<C, S>> {
    let internal = convert_shell_in(shell).ok_or(ShellError::UnsupportedGeometry {
        context: "only planar surfaces and straight edges are supported",
    })?;
    let solid = thicken_internal(&internal, thickness)?;
    Ok(convert_solid_out(&solid))
}

/// Hollows a planar-faced solid to a uniform wall `thickness`, optionally
/// removing the faces in `open_faces`.
///
/// With no open faces the result is a closed hollow solid with two boundary
/// shells (outer and inner cavity). With open faces, each open face becomes
/// a ring (its original boundary plus the offset boundary as a hole) and the
/// cavity connects to the outside through it; open-face planes still
/// constrain their vertices, so the cavity reaches all the way to the open
/// face's plane.
///
/// # Errors
/// - [`ShellError::UnsupportedSolidBoundaries`]: the solid is already hollow
///   (more than one boundary shell).
/// - [`ShellError::FaceNotFound`]: an id in `open_faces` is not a face of
///   the solid.
/// - [`ShellError::AdjacentOpenFaces`]: two open faces share an edge.
/// - plus every error listed for [`thicken`].
pub fn hollow<C: ThickenableCurve, S: ThickenableSurface>(
    solid: &ExternalSolid<C, S>,
    thickness: f64,
    open_faces: &[ExternalFaceId<S>],
) -> Result<ExternalSolid<C, S>> {
    let boundaries = solid.boundaries();
    let [external_shell] = boundaries.as_slice() else {
        return Err(ShellError::UnsupportedSolidBoundaries(boundaries.len()));
    };
    let open_indices = open_faces
        .iter()
        .map(|face_id| {
            external_shell
                .iter()
                .position(|face| face.id() == *face_id)
                .ok_or(ShellError::FaceNotFound)
        })
        .collect::<Result<HashSet<usize>>>()?;
    let internal = convert_shell_in(external_shell).ok_or(ShellError::UnsupportedGeometry {
        context: "only planar surfaces and straight edges are supported",
    })?;
    let solid = hollow_internal(&internal, thickness, &open_indices)?;
    Ok(convert_solid_out(&solid))
}

fn validate_thickness(thickness: f64) -> Result<()> {
    if thickness <= TOLERANCE {
        Err(ShellError::NonPositiveThickness(thickness))
    } else {
        Ok(())
    }
}

fn thicken_internal(shell: &Shell, thickness: f64) -> Result<Solid> {
    validate_thickness(thickness)?;
    if !shell.is_connected() {
        return Err(ShellError::NotOrientable {
            context: "shell is not connected",
        });
    }
    let closed = match shell.shell_condition() {
        ShellCondition::Closed => true,
        ShellCondition::Oriented => false,
        _ => {
            return Err(ShellError::NotOrientable {
                context: "an edge is shared by more than two faces or with incompatible orientation",
            });
        }
    };
    let elements = offset_elements(shell, thickness, &HashSet::default())?;
    let boundaries = if closed {
        let inner: Shell = elements.faces.iter().map(Face::inverse).collect();
        vec![shell.clone(), inner]
    } else {
        let walls = boundary_walls(shell, &elements)?;
        let combined: Shell = shell
            .iter()
            .cloned()
            .chain(elements.faces.iter().map(Face::inverse))
            .chain(walls)
            .collect();
        vec![combined]
    };
    Solid::try_new(boundaries).map_err(|source| ShellError::InvalidOutputTopology { source })
}

fn hollow_internal(shell: &Shell, thickness: f64, open_faces: &HashSet<usize>) -> Result<Solid> {
    if open_faces.is_empty() {
        return thicken_internal(shell, thickness);
    }
    validate_thickness(thickness)?;
    let mut open_edge_ids: HashSet<EdgeId> = HashSet::default();
    for &face_index in open_faces {
        for edge in shell[face_index].edge_iter() {
            if !open_edge_ids.insert(edge.id()) {
                return Err(ShellError::AdjacentOpenFaces);
            }
        }
    }
    let elements = offset_elements(shell, thickness, open_faces)?;
    let combined = shell
        .iter()
        .enumerate()
        .map(
            |(face_index, face)| match open_faces.contains(&face_index) {
                true => Ok(vec![ring_face(face, &elements)?]),
                false => Ok(vec![face.clone(), elements.faces[face_index].inverse()]),
            },
        )
        .collect::<Result<Vec<_>>>()?
        .into_iter()
        .flatten()
        .collect::<Shell>();
    Solid::try_new(vec![combined]).map_err(|source| ShellError::InvalidOutputTopology { source })
}

/// Builds the ring face replacing an open face: the original boundary loops
/// plus the offset boundary loops, inverted, as holes.
fn ring_face(face: &Face, elements: &OffsetElements) -> Result<Face> {
    let holes = face
        .boundaries()
        .iter()
        .map(|wire| {
            let offset_wire = wire
                .edge_iter()
                .map(|edge| elements.oriented_edge(edge))
                .collect::<Result<Wire>>()?;
            Ok(offset_wire.inverse())
        })
        .collect::<Result<Vec<_>>>()?;
    let boundaries = face.boundaries().into_iter().chain(holes).collect();
    Face::try_new(boundaries, face.oriented_surface())
        .map_err(|source| ShellError::InvalidOutputTopology { source })
}

/// Builds planar side walls along the boundary wires of an open shell.
///
/// For each boundary edge `v0 -> v1` (oriented as in its face, interior on
/// the left), the wall is the quad `v1 -> v0 -> v0' -> v1'` on the plane
/// spanned by the edge and the offset direction, so its normal points away
/// from the thickened material.
fn boundary_walls(shell: &Shell, elements: &OffsetElements) -> Result<Vec<Face>> {
    let mut side_edges: HashMap<VertexId, Edge> = HashMap::default();
    let mut side_edge = |vertex: &Vertex| -> Result<Edge> {
        match side_edges.get(&vertex.id()) {
            Some(edge) => Ok(edge.clone()),
            None => {
                let offset_vertex = elements.vertex(vertex)?;
                let line = Line(vertex.point(), offset_vertex.point());
                let edge = Edge::new(vertex, offset_vertex, line);
                side_edges.insert(vertex.id(), edge.clone());
                Ok(edge)
            }
        }
    };
    shell
        .extract_boundaries()
        .iter()
        .flat_map(Wire::edge_iter)
        .map(|edge| {
            let (v0, v1) = edge.ends();
            let offset_edge = elements.oriented_edge(edge)?;
            let (side0, side1) = (side_edge(v0)?, side_edge(v1)?);
            let (p0, p1) = (v0.point(), v1.point());
            let (q0, q1) = (offset_edge.front().point(), offset_edge.back().point());
            let normal = (p0 - p1).cross(q0 - p1);
            if normal.magnitude() < TOLERANCE {
                return Err(ShellError::GeometryFailed {
                    context: "degenerate side wall",
                });
            }
            if (q1 - p1).dot(normal.normalize()).abs() > TOLERANCE {
                return Err(ShellError::UnsupportedGeometry {
                    context: "boundary vertices offset out of the side-wall plane",
                });
            }
            let wall_wire: Wire = [edge.inverse(), side0, offset_edge, side1.inverse()]
                .into_iter()
                .collect();
            Face::try_new(vec![wall_wire], Plane::new(p1, p0, q0))
                .map_err(|source| ShellError::InvalidOutputTopology { source })
        })
        .collect()
}
