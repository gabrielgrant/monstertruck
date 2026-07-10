use monstertruck_geometry::prelude::*;
use rustc_hash::FxHashSet as HashSet;

use super::error::DraftError;
use super::offset::drafted_shell;
use super::types::*;

type Result<T> = std::result::Result<T, DraftError>;
type ExternalSolid<C, S> = monstertruck_topology::Solid<Point3, C, S>;
type ExternalFaceId<S> = monstertruck_topology::FaceId<S>;

/// Draft angles are restricted to the open interval `(0°, 90°)`.
const MAX_ANGLE: Rad<f64> = Rad(std::f64::consts::FRAC_PI_2);

/// Applies a planar draft (mold taper) to the selected `faces` of `solid`.
///
/// Each selected face is a plane; it is rotated by `angle` about its
/// intersection line with `neutral_plane`, tilting so the face opens toward
/// `pull_direction` as it moves away from the neutral plane -- the standard
/// injection-molding draft, which eases release of the part from its mold
/// along the pull direction. Faces not in `faces` keep their original
/// planes. Every vertex and edge of the solid's boundary shell is then
/// reconciled by intersecting the (possibly unchanged) planes of its
/// adjacent faces, exactly as in the shell/thicken module's offset
/// reconciliation: draft's per-face motion is a rotation about a line
/// instead of a translation along a normal, but the shared-vertex problem it
/// produces is the same plane-intersection problem, solved the same way.
///
/// # Errors
/// - [`DraftError::UnsupportedGeometry`]: a surface is not a
///   [`Plane`](monstertruck_geometry::prelude::Plane), or an edge curve is
///   not a [`Line`](monstertruck_geometry::prelude::Line).
/// - [`DraftError::UnsupportedSolidBoundaries`]: the solid is hollow (more
///   than one boundary shell).
/// - [`DraftError::FaceNotFound`]: an id in `faces` is not a face of the
///   solid.
/// - [`DraftError::NoFacesSelected`]: `faces` is empty.
/// - [`DraftError::InvalidAngle`]: `angle` is not strictly between 0 and 90
///   degrees.
/// - [`DraftError::DegeneratePullDirection`]: `pull_direction` is the zero
///   vector.
/// - [`DraftError::ParallelToNeutralPlane`]: a selected face is parallel to
///   `neutral_plane` and so has no intersection line to rotate about.
/// - [`DraftError::AmbiguousPullDirection`]: `pull_direction` is too close
///   to parallel with a selected face's hinge line to determine which way
///   "toward the pull direction" tilts that face.
/// - [`DraftError::ExcessiveAngle`]: `angle` collapses or inverts an edge.
/// - [`DraftError::DegenerateVertex`] / [`DraftError::UnsupportedVertexDegree`]:
///   a vertex cannot be reconciled by plane intersection.
pub fn draft<C: DraftableCurve, S: DraftableSurface, A: Into<Rad<f64>>>(
    solid: &ExternalSolid<C, S>,
    faces: &[ExternalFaceId<S>],
    neutral_plane: Plane,
    pull_direction: Vector3,
    angle: A,
) -> Result<ExternalSolid<C, S>> {
    let angle = angle.into();
    if angle.0 <= 0.0 || angle.0 >= MAX_ANGLE.0 {
        return Err(DraftError::InvalidAngle { radians: angle.0 });
    }
    if pull_direction.magnitude() < TOLERANCE {
        return Err(DraftError::DegeneratePullDirection);
    }
    if faces.is_empty() {
        return Err(DraftError::NoFacesSelected);
    }
    let boundaries = solid.boundaries();
    let [external_shell] = boundaries.as_slice() else {
        return Err(DraftError::UnsupportedSolidBoundaries(boundaries.len()));
    };
    let drafted_indices = faces
        .iter()
        .map(|face_id| {
            external_shell
                .iter()
                .position(|face| face.id() == *face_id)
                .ok_or(DraftError::FaceNotFound)
        })
        .collect::<Result<HashSet<usize>>>()?;
    let internal = convert_shell_in(external_shell).ok_or(DraftError::UnsupportedGeometry {
        context: "only planar surfaces and straight edges are supported",
    })?;
    let shell = drafted_shell(
        &internal,
        &drafted_indices,
        neutral_plane,
        pull_direction,
        angle,
    )?;
    let solid = Solid::try_new(vec![shell])
        .map_err(|source| DraftError::InvalidOutputTopology { source })?;
    Ok(convert_solid_out(&solid))
}
