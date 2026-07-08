//! Fillet operations for [`Shell`](monstertruck_topology::Shell) edges.
//!
//! Provides rolling-ball fillet operations: single-edge fillets,
//! fillets with side face updates, and fillets along open or closed wire chains.
//! The [`fillet_edges`] function provides a high-level API for selected edges.
//! The [`fillet_edges_by_id`] function resolves face adjacency from edge IDs.
//!
//! ## Batch-filleting edges that share a vertex
//!
//! [`fillet_edges`]/[`fillet_edges_by_id`] group a multi-edge selection into
//! chains of contiguous, vertex-sharing edges on a single face boundary
//! (`group_edges_into_chains`), then fillet each chain in turn as one wire
//! via [`fillet_along_wire`], re-resolving edge identities between chains
//! since an earlier chain's mutation can split or replace edges an
//! as-yet-unprocessed chain still refers to by their original ID
//! (`rematch_selected_edge_id`). If a chain's edges fail to re-resolve down
//! to fewer than 2 edges, `fillet_along_wire` cannot find a shared face for a
//! degenerate wire and falls back to filleting the chain's edges one at a
//! time (`fillet_edges_by_id`'s `wire_fillet_failed` branch); a chain that
//! genuinely cannot be resolved is rolled back to its pre-chain state rather
//! than left partially filleted. This keeps ordinary "round these two edges"
//! selections -- including two edges that meet at a corner -- both panic-free
//! and manifold.
//!
//! `fillet_along_wire`'s closed- and open-wire face construction
//! (`fillet_along_wire_closed`/`_open`) has a separate, pre-existing
//! limitation: for a wire chain that spans three or more consecutive side
//! faces (e.g. a closed loop around all four edges of a cuboid's top face),
//! each side face's shared vertical edge is cut independently by
//! `cut_face_by_last_bezier`/`cut_face_by_bezier` once per adjacent side
//! face, rather than being reconciled into one edge shared by both faces.
//! The two independently-cut copies are geometrically coincident but are
//! distinct [`Edge`](monstertruck_topology::Edge) objects, so the result is
//! geometrically correct but not strictly manifold
//! (`Shell::shell_condition` reports `Oriented`, not `Closed`, at those
//! seams) -- see `fillet_edges_cuboid_top_4`, `fillet_edges_cuboid_top_and_bottom`,
//! and `fillet_edges_all_twelve` in the test module. Closing that gap needs
//! the two side faces that share such a vertical edge to agree on a single
//! cut result (a small cross-face reconciliation pass, or genuine
//! corner-blend surfaces at 3+-edge junctions); neither exists yet, so
//! batch-filleting an entire closed loop (as opposed to a couple of edges at
//! one corner) can still leave zero-width, geometrically-coincident seams
//! rather than a single shared edge. `fillet_edges`/`fillet_edges_by_id`
//! never panic and never produce an `Irregular` shell (an edge shared by
//! more than two faces) for this case, but they do not yet guarantee
//! `Closed`.

#[allow(private_interfaces)]
mod edge_select;
#[allow(private_interfaces)]
mod ops;

mod convert;
mod error;
mod geometry;
mod params;
mod topology;
mod types;

#[cfg(test)]
mod tests;

pub use convert::{FilletIntersectionCurve, FilletableCurve, FilletableSurface};
pub use edge_select::{fillet_edges, fillet_edges_by_id};
pub use error::FilletError;
pub use ops::{fillet, fillet_along_wire, fillet_with_side};
pub use params::{FilletOptions, FilletProfile, FilletRadius};
pub use types::ParameterCurveLinear;
