use monstertruck_topology::errors::Error as TopologyError;
use thiserror::Error;

use crate::planar::PlanarError;

/// Errors that can occur during shell/thicken operations.
#[derive(Debug, Error)]
pub enum ShellError {
    /// The shell contains geometry that cannot be handled by the planar-only
    /// implementation: a surface that is not a [`Plane`](monstertruck_geometry::prelude::Plane)
    /// or an edge whose curve is not a [`Line`](monstertruck_geometry::prelude::Line).
    #[error("Unsupported geometry: {context}.")]
    UnsupportedGeometry {
        /// Description of which piece of geometry could not be handled.
        context: &'static str,
    },
    /// The shell is not a single connected, orientable surface (it is
    /// disconnected, non-manifold, or its faces disagree on orientation
    /// across a shared edge).
    #[error("Shell is not manifold and orientable: {context}.")]
    NotOrientable {
        /// Description of the offending condition.
        context: &'static str,
    },
    /// One of the requested open faces was not found in the solid.
    #[error("Face not found in solid.")]
    FaceNotFound,
    /// [`hollow`](super::hollow) requires a solid with exactly one boundary
    /// shell (a solid that is not already hollow).
    #[error("Solid must have exactly one boundary shell, found {0}.")]
    UnsupportedSolidBoundaries(usize),
    /// The requested thickness is not positive.
    #[error("Thickness must be positive, got {0}.")]
    NonPositiveThickness(f64),
    /// The requested thickness exceeds the local feature size: offsetting
    /// collapses or inverts at least one edge. This local check does not
    /// detect non-local self-intersections -- see the
    /// [module documentation](super) for details.
    #[error("Thickness {thickness} is too large: offsetting collapses or inverts an edge.")]
    ThicknessTooLarge {
        /// The requested thickness.
        thickness: f64,
    },
    /// A vertex offset could not be computed because the adjacent offset
    /// planes do not intersect in a single well-defined point (nearly
    /// parallel or anti-parallel planes, or coplanar faces with conflicting
    /// offsets).
    #[error("Degenerate vertex: adjacent offset planes do not intersect cleanly.")]
    DegenerateVertex,
    /// A vertex is adjacent to more distinct face planes than this
    /// constrained implementation supports (more than 3).
    #[error("Vertex adjacent to {0} distinct face planes; only up to 3 is supported.")]
    UnsupportedVertexDegree(usize),
    /// Two of the requested open faces share an edge. Offsetting collapses
    /// the shared edge onto itself, so adjacent open faces are rejected.
    #[error("Two open faces share an edge; adjacent open faces are not supported.")]
    AdjacentOpenFaces,
    /// A geometry computation failed unexpectedly.
    #[error("Geometry failed: {context}.")]
    GeometryFailed {
        /// Description of which step failed.
        context: &'static str,
    },
    /// The generated shell or solid is topologically invalid.
    #[error("Invalid output topology: {source}.")]
    InvalidOutputTopology {
        /// Topology validation error.
        #[source]
        source: TopologyError,
    },
}

impl From<PlanarError> for ShellError {
    fn from(error: PlanarError) -> Self {
        match error {
            PlanarError::DegenerateVertex => ShellError::DegenerateVertex,
            PlanarError::UnsupportedVertexDegree(degree) => {
                ShellError::UnsupportedVertexDegree(degree)
            }
        }
    }
}
