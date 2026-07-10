use monstertruck_topology::errors::Error as TopologyError;
use thiserror::Error;

/// Errors that can occur during draft (mold-taper) operations.
#[derive(Debug, Error)]
pub enum DraftError {
    /// The shell contains geometry that cannot be handled by the
    /// planar-only implementation: a surface that is not a
    /// [`Plane`](monstertruck_geometry::prelude::Plane) or an edge whose
    /// curve is not a [`Line`](monstertruck_geometry::prelude::Line).
    #[error("Unsupported geometry: {context}.")]
    UnsupportedGeometry {
        /// Description of which piece of geometry could not be handled.
        context: &'static str,
    },
    /// One of the requested faces was not found in the solid.
    #[error("Face not found in solid.")]
    FaceNotFound,
    /// [`draft`](super::draft) requires a solid with exactly one boundary
    /// shell (a solid that is not hollow).
    #[error("Solid must have exactly one boundary shell, found {0}.")]
    UnsupportedSolidBoundaries(usize),
    /// No faces were selected to draft.
    #[error("No faces were selected for draft.")]
    NoFacesSelected,
    /// The requested angle is not strictly between 0 and 90 degrees.
    #[error("Draft angle must be strictly between 0 and 90 degrees, got {radians} rad.")]
    InvalidAngle {
        /// The requested angle, in radians.
        radians: f64,
    },
    /// The pull direction was the zero vector.
    #[error("Pull direction must be a nonzero vector.")]
    DegeneratePullDirection,
    /// A selected face's plane is parallel to the neutral plane, so it has
    /// no intersection line to hinge about.
    #[error("Face is parallel to the neutral plane: it has no intersection line to draft about.")]
    ParallelToNeutralPlane,
    /// The pull direction is too close to parallel with a selected face's
    /// hinge line to determine which way the face should tilt.
    #[error("Pull direction is ambiguous for a face's hinge line (nearly parallel to it).")]
    AmbiguousPullDirection,
    /// The requested angle exceeds the local feature size: drafting
    /// collapses or inverts at least one edge. This local check does not
    /// detect non-local self-intersections -- see the
    /// [module documentation](super) for details.
    #[error("Draft angle {radians} rad is too large: it collapses or inverts an edge.")]
    ExcessiveAngle {
        /// The requested angle, in radians.
        radians: f64,
    },
    /// A vertex reconciliation could not be computed because the adjacent
    /// drafted planes do not intersect in a single well-defined point
    /// (nearly parallel or anti-parallel planes, or coplanar faces with
    /// conflicting planes).
    #[error("Degenerate vertex: adjacent drafted planes do not intersect cleanly.")]
    DegenerateVertex,
    /// A vertex is adjacent to more distinct face planes than this
    /// constrained implementation supports (more than 3).
    #[error("Vertex adjacent to {0} distinct face planes; only up to 3 is supported.")]
    UnsupportedVertexDegree(usize),
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
