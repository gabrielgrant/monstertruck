//! Planar draft (mold taper) about a neutral plane.
//!
//! [`draft`] tilts selected planar faces of a solid by a fixed angle about
//! their intersection line with a neutral plane, the way an injection mold
//! tapers a part's walls so it releases cleanly along the pull direction.
//! Each selected face's plane is replaced by that plane rotated about
//! `face plane ∩ neutral plane`; faces not selected keep their original
//! planes. Shared vertices and edges are then rebuilt by intersecting the
//! (possibly unchanged) planes of every adjacent face -- the same
//! plane-reconciliation approach the shell/thicken module uses for
//! offsetting (see that module's docs), just with a rotated plane per face
//! in place of a translated one. Every original vertex and edge maps to
//! exactly one drafted vertex/edge, so rebuilt faces share their edges by
//! construction and results are validated strictly manifold by
//! [`Solid::try_new`](monstertruck_topology::Solid::try_new).
//!
//! Unlike shell/thicken, draft never adds or removes faces: it replaces the
//! solid's single boundary shell with a shell of the same face count and
//! connectivity, only vertices and edges move.
//!
//! ## Supported geometry
//!
//! Deliberately as constrained as the shell/thicken module's first version,
//! for the same reasons:
//!
//! - **Planar faces and straight edges only.** Surfaces must convert to
//!   [`Plane`](monstertruck_geometry::prelude::Plane) and edge curves to
//!   [`Line`](monstertruck_geometry::prelude::Line) (see [`DraftableSurface`]
//!   and [`DraftableCurve`]), checked by type, not by sampling.
//! - **Vertices adjacent to at most 3 distinct face planes.** A vertex where
//!   4 or more planes meet has no single reconciled point in general and is
//!   rejected with [`DraftError::UnsupportedVertexDegree`].
//! - **Selected faces must not be parallel to the neutral plane**
//!   ([`DraftError::ParallelToNeutralPlane`]): a face parallel to the
//!   neutral plane has no intersection line to hinge about.
//! - **The pull direction must disambiguate the tilt.** If it is too close
//!   to parallel with a selected face's hinge line, which way "toward the
//!   pull direction" tilts that face is undefined
//!   ([`DraftError::AmbiguousPullDirection`]).
//! - **The draft angle must invert nothing.** An angle that collapses or
//!   inverts an edge is rejected ([`DraftError::ExcessiveAngle`]) rather
//!   than producing a self-intersecting shell; angles are further
//!   restricted to the open interval `(0°, 90°)`
//!   ([`DraftError::InvalidAngle`]).
//! - **Local self-intersection check only**, as in shell/thicken: two
//!   drafted walls far apart in the topology crossing in space is not
//!   detected; the caller must keep the angle and part geometry sane.

mod error;
mod offset;
mod ops;
mod types;

#[cfg(test)]
mod tests;

pub use error::DraftError;
pub use ops::draft;
pub use types::{DraftableCurve, DraftableSurface};
