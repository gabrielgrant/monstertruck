//! Crate for operation shapes. Provides boolean operations to Solid, and shape healing for importing shapes from other CAD systems.

#![cfg_attr(not(debug_assertions), deny(warnings))]
#![deny(clippy::all, rust_2018_idioms)]
#![warn(
    missing_docs,
    missing_debug_implementations,
    trivial_casts,
    trivial_numeric_casts,
    unsafe_code,
    unstable_features,
    unused_import_braces,
    unused_qualifications
)]

mod healing;
pub use healing::{RobustSplitClosedEdgesAndFaces, SplitClosedEdgesAndFaces, extract_healed};
mod transversal;
pub use transversal::{
    ShapeOpsCurve, ShapeOpsError, ShapeOpsSurface, and, difference, or, symmetric_difference,
};
mod alternative;
mod planar;
pub mod shell;
pub use shell::{ShellError, ThickenableCurve, ThickenableSurface, hollow, thicken};
pub mod fillet;
pub use fillet::{
    FilletError, FilletIntersectionCurve, FilletOptions, FilletProfile, FilletRadius,
    FilletableCurve, FilletableSurface, ParameterCurveLinear, fillet, fillet_along_wire,
    fillet_edges, fillet_edges_by_id, fillet_with_side,
};
pub mod draft;
pub use draft::{DraftError, DraftableCurve, DraftableSurface, draft};
