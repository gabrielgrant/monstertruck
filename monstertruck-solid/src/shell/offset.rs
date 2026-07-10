use monstertruck_geometry::prelude::*;
use rustc_hash::{FxHashMap as HashMap, FxHashSet as HashSet};
use smallvec::SmallVec;

use super::error::ShellError;
use super::types::*;

type Result<T> = std::result::Result<T, ShellError>;

/// A plane constraint on a vertex offset: the offset displacement `w` must
/// satisfy `normal . w = -offset` (move `offset` against the face normal).
type Constraint = (Vector3, f64);

/// Offset images of every vertex, edge, and face of a planar shell.
///
/// Every original vertex maps to exactly one offset vertex and every
/// original edge to exactly one offset edge, so faces rebuilt from these
/// maps share their edges by construction (single [`Edge`] objects, no
/// coincident duplicates).
pub(super) struct OffsetElements {
    /// Offset vertex for each original vertex.
    pub(super) vertices: HashMap<VertexId, Vertex>,
    /// Offset edge for each original edge, in absolute orientation.
    pub(super) edges: HashMap<EdgeId, Edge>,
    /// Offset face for each original face (same index and orientation).
    pub(super) faces: Vec<Face>,
}

impl OffsetElements {
    /// Returns the offset vertex for `vertex`.
    pub(super) fn vertex(&self, vertex: &Vertex) -> Result<&Vertex> {
        self.vertices
            .get(&vertex.id())
            .ok_or(ShellError::GeometryFailed {
                context: "offset vertex lookup",
            })
    }

    /// Returns the offset edge for `edge`, oriented like `edge`.
    pub(super) fn oriented_edge(&self, edge: &Edge) -> Result<Edge> {
        let offset_edge = self
            .edges
            .get(&edge.id())
            .ok_or(ShellError::GeometryFailed {
                context: "offset edge lookup",
            })?;
        Ok(match edge.orientation() {
            true => offset_edge.clone(),
            false => offset_edge.inverse(),
        })
    }
}

/// Solves for the displacement `w` with `n_i . w = -t_i` for every
/// constraint `(n_i, t_i)`.
///
/// With one constraint the unique minimal solution `-t n` is used; with two,
/// the solution in `span(n_0, n_1)` (so `w` stays perpendicular to the
/// shared edge of the two faces); with three, the unique intersection of the
/// three offset planes.
fn solve_offset(constraints: &[Constraint]) -> Result<Vector3> {
    match *constraints {
        [(normal, offset)] => Ok(-offset * normal),
        [(normal0, offset0), (normal1, offset1)] => {
            let cos = normal0.dot(normal1);
            let det = 1.0 - cos * cos;
            if det.abs() < TOLERANCE2 {
                Err(ShellError::DegenerateVertex)
            } else {
                let a = (cos * offset1 - offset0) / det;
                let b = (cos * offset0 - offset1) / det;
                Ok(a * normal0 + b * normal1)
            }
        }
        [(normal0, offset0), (normal1, offset1), (normal2, offset2)] => {
            let matrix = Matrix3::from_cols(normal0, normal1, normal2).transpose();
            if matrix.determinant().abs() < TOLERANCE {
                Err(ShellError::DegenerateVertex)
            } else {
                // SAFETY: the determinant was just checked to be non-zero.
                let inverse = matrix.invert().unwrap();
                Ok(inverse * -Vector3::new(offset0, offset1, offset2))
            }
        }
        _ => Err(ShellError::UnsupportedVertexDegree(constraints.len())),
    }
}

/// Collects the distinct plane constraints of a vertex, merging near-equal
/// normals and rejecting coplanar faces that request conflicting offsets.
fn vertex_constraints(
    face_data: &[Constraint],
    face_indices: &[usize],
) -> Result<SmallVec<[Constraint; 4]>> {
    let mut constraints: SmallVec<[Constraint; 4]> = SmallVec::new();
    for &face_index in face_indices {
        let (normal, offset) = face_data[face_index];
        match constraints
            .iter()
            .find(|&&(other, _)| (other - normal).magnitude2() < TOLERANCE2)
        {
            Some(&(_, other_offset)) if (other_offset - offset).abs() < TOLERANCE => {}
            Some(_) => return Err(ShellError::DegenerateVertex),
            None => constraints.push((normal, offset)),
        }
    }
    Ok(constraints)
}

/// Computes the offset image of every vertex, edge, and face of `shell`.
///
/// Each face is offset by `thickness` against its oriented normal, except
/// faces whose index is in `open_faces`, which are offset by zero (they
/// still constrain their vertices to slide within their own plane). Shared
/// vertices are reconciled by intersecting the adjacent offset planes.
pub(super) fn offset_elements(
    shell: &Shell,
    thickness: f64,
    open_faces: &HashSet<usize>,
) -> Result<OffsetElements> {
    let face_data: Vec<Constraint> = shell
        .iter()
        .enumerate()
        .map(|(face_index, face)| {
            let offset = match open_faces.contains(&face_index) {
                true => 0.0,
                false => thickness,
            };
            (face.oriented_surface().normal(), offset)
        })
        .collect();

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
            let constraints = vertex_constraints(&face_data, &face_indices)?;
            let displacement = solve_offset(&constraints)?;
            Ok((vertex_id, Vertex::new(vertex.point() + displacement)))
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
                .ok_or(ShellError::GeometryFailed {
                    context: "offset vertex lookup for edge",
                })
        };
        let (new_front, new_back) = (lookup(front)?, lookup(back)?);
        let original = back.point() - front.point();
        let offset = new_back.point() - new_front.point();
        if offset.magnitude() < TOLERANCE || offset.dot(original) <= 0.0 {
            return Err(ShellError::ThicknessTooLarge { thickness });
        }
        let line = Line(new_front.point(), new_back.point());
        edges.insert(edge.id(), Edge::new(new_front, new_back, line));
    }

    let faces = shell
        .iter()
        .enumerate()
        .map(|(face_index, face)| {
            let (normal, offset) = face_data[face_index];
            let displacement = -offset * normal;
            let boundaries = face
                .absolute_boundaries()
                .iter()
                .map(|wire| {
                    wire.iter()
                        .map(|edge| {
                            let new_edge =
                                edges.get(&edge.id()).ok_or(ShellError::GeometryFailed {
                                    context: "offset edge lookup for face",
                                })?;
                            Ok(match edge.orientation() {
                                true => new_edge.clone(),
                                false => new_edge.inverse(),
                            })
                        })
                        .collect::<Result<Wire>>()
                })
                .collect::<Result<Vec<_>>>()?;
            let surface = face.surface();
            let origin = surface.origin() + displacement;
            let plane = Plane::new(origin, origin + surface.u_axis(), origin + surface.v_axis());
            let mut new_face = Face::try_new(boundaries, plane)
                .map_err(|source| ShellError::InvalidOutputTopology { source })?;
            if !face.orientation() {
                new_face.invert();
            }
            Ok(new_face)
        })
        .collect::<Result<Vec<_>>>()?;

    Ok(OffsetElements {
        vertices,
        edges,
        faces,
    })
}
