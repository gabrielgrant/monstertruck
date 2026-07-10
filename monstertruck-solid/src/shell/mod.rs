//! Shell/thicken operations for planar-faced shells and solids.
//!
//! [`thicken`] turns an open or closed [`Shell`](monstertruck_topology::Shell)
//! into a solid of uniform wall thickness; [`hollow`] shells out a solid,
//! optionally opening selected faces. Both offset every face plane by the
//! wall thickness against its oriented normal and rebuild shared edges and
//! vertices by intersecting the adjacent offset planes -- convex and concave
//! corners alike. Every original vertex and edge maps to exactly one offset
//! vertex/edge, so the offset faces share single [`Edge`](monstertruck_topology::Edge)
//! objects by construction and results are strictly manifold
//! (`Shell::shell_condition` reports `Closed`; [`Solid::try_new`](monstertruck_topology::Solid::try_new)
//! validates this on every result). This differs from the fillet module's
//! per-face independent cuts, which can leave geometrically-coincident
//! duplicate seam edges (see the fillet module docs); the plane-intersection
//! reconciliation here is only possible because the geometry is restricted
//! to planes.
//!
//! ## Supported geometry
//!
//! This is a deliberately constrained first version:
//!
//! - **Planar faces and straight edges only.** Surfaces must convert to
//!   [`Plane`](monstertruck_geometry::prelude::Plane) and edge curves to
//!   [`Line`](monstertruck_geometry::prelude::Line) (see [`ThickenableSurface`]
//!   and [`ThickenableCurve`]). The check is by type, not by sampling: a
//!   NURBS surface that happens to be flat is still rejected with
//!   [`ShellError::UnsupportedGeometry`]. Curved faces would need
//!   `OffsetSurface` geometry plus surface-surface intersection for the
//!   corners; neither is wired up here.
//! - **Vertices adjacent to at most 3 distinct face planes.** A vertex where
//!   4 or more planes meet has no single offset point in general and is
//!   rejected with [`ShellError::UnsupportedVertexDegree`]. Coplanar
//!   neighbors are merged before counting.
//! - **Local self-intersection check only.** A thickness that collapses or
//!   inverts an edge (e.g. more than half the width of a box) is rejected
//!   with [`ShellError::ThicknessTooLarge`]. Non-local self-intersection --
//!   two walls far apart in the topology crossing in space -- is **not**
//!   detected; the caller must keep the thickness below the model's minimal
//!   feature size.
//! - **Open faces of [`hollow`] must not share edges** with each other
//!   ([`ShellError::AdjacentOpenFaces`]): the shared edge would offset onto
//!   its own line and produce a degenerate ring face.
//! - **Side walls of an open-shell [`thicken`] must be planar.** Boundary
//!   vertices offset within the span of their adjacent face normals; if the
//!   two ends of a boundary edge offset out of a common plane (possible on
//!   non-flat open shells with oblique corners), the wall cannot be a single
//!   planar face and the operation fails with
//!   [`ShellError::UnsupportedGeometry`] rather than producing a warped
//!   quad.

mod error;
mod offset;
mod ops;
mod types;

#[cfg(test)]
mod tests;

pub use error::ShellError;
pub use ops::{hollow, thicken};
pub use types::{ThickenableCurve, ThickenableSurface};
