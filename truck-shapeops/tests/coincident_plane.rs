//! MCVE mirroring openshape's fs-kernel-truck regression matrix
//! (crates/fs-kernel-truck/tests/boolean_union_matrix.rs and
//! upstream_tol_sweep.rs in https://github.com/.../openshape) directly
//! against `truck_shapeops::and`/`or`, bypassing openshape entirely, using
//! only `truck_modeling` builder cuboids.
//!
//! Confirms: touching-face / coincident-plane cuboid pairs make `or`/`and`
//! return `None` at every tested tolerance, and a nested-corner cuboid pair
//! (three coincident planes meeting at a shared corner) panics inside
//! truck-topology with "This wire is not simple" rather than returning
//! `None`.

use truck_modeling::{builder, Point3, Solid, Vector3};

fn cuboid(lo: [f64; 3], hi: [f64; 3]) -> Solid {
    let v = builder::vertex(Point3::new(lo[0], lo[1], lo[2]));
    let e = builder::tsweep(&v, Vector3::new(hi[0] - lo[0], 0.0, 0.0));
    let f = builder::tsweep(&e, Vector3::new(0.0, hi[1] - lo[1], 0.0));
    builder::tsweep(&f, Vector3::new(0.0, 0.0, hi[2] - lo[2]))
}

const TOLS: [f64; 7] = [1e-6, 1e-4, 1e-3, 0.01, 0.05, 0.1, 0.5];

/// (a) Clean overlap, no coincident faces: `or` should succeed at every
/// tolerance. Sanity check / control case.
#[test]
fn a_clean_overlap_succeeds() {
    let a = cuboid([0.0, 0.0, 0.0], [4.0, 4.0, 4.0]);
    let b = cuboid([2.0, 2.0, 2.0], [6.0, 6.0, 6.0]);
    for tol in TOLS {
        let result = truck_shapeops::or(&a, &b, tol);
        assert!(result.is_some(), "clean overlap should succeed at tol={tol}");
    }
}

/// (b) Touching faces, cubes share the x=2 plane exactly (zero-volume
/// overlap). CONFIRMED: `or` returns `None` at every tolerance from 1e-6 to
/// 0.5.
#[test]
fn b_touching_face_fails_at_every_tolerance() {
    let a = cuboid([0.0, 0.0, 0.0], [2.0, 2.0, 2.0]);
    let b = cuboid([2.0, 0.0, 0.0], [4.0, 2.0, 2.0]);
    for tol in TOLS {
        let result = truck_shapeops::or(&a, &b, tol);
        assert!(
            result.is_none(),
            "touching-face union expected None (documented degeneracy), got Some at tol={tol}"
        );
    }
}

/// Same geometry as (b), but for `and` -- confirms it's not `or`-specific.
#[test]
fn b_touching_face_and_also_fails() {
    let a = cuboid([0.0, 0.0, 0.0], [2.0, 2.0, 2.0]);
    let b = cuboid([2.0, 0.0, 0.0], [4.0, 2.0, 2.0]);
    for tol in TOLS {
        let result = truck_shapeops::and(&a, &b, tol);
        assert!(
            result.is_none(),
            "touching-face intersection expected None, got Some at tol={tol}"
        );
    }
}

/// (i) Tiny gap (0.0001): genuinely disjoint at that scale -- works (or
/// returns Some) since it no longer shares a coincident plane.
#[test]
fn i_near_touching_epsilon_gap_succeeds() {
    let a = cuboid([0.0, 0.0, 0.0], [2.0, 2.0, 2.0]);
    let b = cuboid([2.0001, 0.0, 0.0], [4.0001, 2.0, 2.0]);
    let result = truck_shapeops::or(&a, &b, 0.05);
    assert!(result.is_some(), "epsilon-gap disjoint cuboids should succeed");
}

/// (j) Tiny overlap (0.0001) instead of exact touching: still fails --
/// confirms this isn't literal floating-point equality, it's a
/// tolerance-scale near-coincidence.
#[test]
fn j_near_touching_epsilon_overlap_fails() {
    let a = cuboid([0.0, 0.0, 0.0], [2.0, 2.0, 2.0]);
    let b = cuboid([1.9999, 0.0, 0.0], [3.9999, 2.0, 2.0]);
    let result = truck_shapeops::or(&a, &b, 0.05);
    assert!(
        result.is_none(),
        "near-touching (0.0001 overlap) union expected None, got Some"
    );
}

/// (f) Nested cubes sharing 3 faces at the same corner (origin): full
/// containment plus 3 coincident planes.
///
/// Before the `divide_one_face` fix (truck-shapeops/src/transversal/divide_face/mod.rs),
/// this panicked inside truck-topology with "This wire is not simple" in
/// debug builds -- and in *release* builds took the `Face::new_unchecked`
/// path instead (via `Face::debug_new`), which performs no validation at
/// all and would have silently produced a face with a self-intersecting
/// boundary rather than erroring. `divide_one_face` now uses
/// `Face::try_new` and propagates failure as `None`, so this is a clean,
/// documented failure (same `Expect::Err`-style outcome as the other
/// coincident-plane cases) in both build profiles instead of a panic or
/// silent corruption.
#[test]
fn f_nested_shared_corner_no_longer_panics() {
    let a = cuboid([0.0, 0.0, 0.0], [2.0, 2.0, 2.0]);
    let b = cuboid([0.0, 0.0, 0.0], [4.0, 4.0, 4.0]);
    let result =
        std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| truck_shapeops::or(&a, &b, 0.05)));
    let result = result.unwrap_or_else(|_| {
        panic!("nested-corner union panicked (regression: divide_one_face's Face::try_new fix should turn this into a clean None)")
    });
    assert!(
        result.is_none(),
        "nested-corner union expected None (documented degeneracy), got Some \
         -- if truck-shapeops fixed coincident-plane booleans upstream, update this test"
    );
}
