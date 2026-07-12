//! Shared plane-intersection solver for planar offset operations.
//!
//! [`shell::offset`](super::shell) and [`draft::offset`](super::draft) both
//! reconcile a shared vertex against the (up to three) distinct face planes
//! it is adjacent to, after each face's plane has separately been replaced --
//! by [`shell`](super::shell)'s per-face scalar offset along the face
//! normal, or [`draft`](super::draft)'s per-face rotation about the neutral
//! plane. Both problems reduce to the same linear system once each face's
//! replacement plane is written as an absolute equation
//! `normal . x = distance`: this module holds that constraint bookkeeping
//! and the 1/2/3-plane solve, generalized from shell/thicken's original
//! "offset every plane by a scalar along its own normal" to draft's
//! "replace each plane with an arbitrary new one" (shell's scalar offsets
//! are just a special case, expressed as absolute equations by
//! `distance = original_distance - offset`).

use monstertruck_geometry::prelude::*;
use smallvec::SmallVec;

/// A plane equation `normal . x = distance`, in absolute (world)
/// coordinates. `normal` is always a unit vector.
pub(crate) type PlaneEq = (Vector3, f64);

/// A vertex's plane-reconciliation step failed.
#[derive(Debug, Clone, Copy)]
pub(crate) enum PlanarError {
    /// The adjacent planes do not intersect in a single well-defined point
    /// (nearly parallel or anti-parallel planes, or coplanar faces that
    /// request conflicting planes).
    DegenerateVertex,
    /// A vertex is adjacent to more distinct face planes than the solver
    /// supports (more than 3).
    UnsupportedVertexDegree(usize),
}

type Result<T> = std::result::Result<T, PlanarError>;

/// Collects the distinct plane constraints of a vertex, merging near-equal
/// normals and rejecting coplanar faces that request conflicting planes.
pub(crate) fn dedup_constraints(
    raw: impl IntoIterator<Item = PlaneEq>,
) -> Result<SmallVec<[PlaneEq; 4]>> {
    let mut constraints: SmallVec<[PlaneEq; 4]> = SmallVec::new();
    for (normal, distance) in raw {
        match constraints
            .iter()
            .find(|&&(other, _)| (other - normal).magnitude2() < TOLERANCE2)
        {
            Some(&(_, other_distance)) if (other_distance - distance).abs() < TOLERANCE => {}
            Some(_) => return Err(PlanarError::DegenerateVertex),
            None => constraints.push((normal, distance)),
        }
    }
    Ok(constraints)
}

/// Finds the point satisfying every plane equation in `constraints`, using
/// `reference` to resolve the under-determined cases (1 or 2 constraints).
///
/// - **1 constraint**: the orthogonal projection of `reference` onto the
///   plane.
/// - **2 constraints**: the point reached from `reference` by a displacement
///   confined to `span(normal0, normal1)` -- perpendicular to the planes'
///   shared line, so a `reference` that lies on both planes' original shared
///   edge keeps that edge's direction.
/// - **3 constraints**: the unique intersection point of the three planes
///   (independent of `reference`).
pub(crate) fn intersect_planes(reference: Point3, constraints: &[PlaneEq]) -> Result<Point3> {
    match *constraints {
        [(normal, distance)] => {
            let excess = normal.dot(reference.to_vec()) - distance;
            Ok(reference - excess * normal)
        }
        [(normal0, distance0), (normal1, distance1)] => {
            let cos = normal0.dot(normal1);
            let det = 1.0 - cos * cos;
            if det.abs() < TOLERANCE2 {
                Err(PlanarError::DegenerateVertex)
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
                Err(PlanarError::DegenerateVertex)
            } else {
                // SAFETY: the determinant was just checked to be non-zero.
                let inverse = matrix.invert().unwrap();
                Ok(Point3::from_vec(
                    inverse * Vector3::new(distance0, distance1, distance2),
                ))
            }
        }
        _ => Err(PlanarError::UnsupportedVertexDegree(constraints.len())),
    }
}
