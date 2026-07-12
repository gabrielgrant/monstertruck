use monstertruck_geometry::prelude::*;
use monstertruck_meshing::prelude::*;
use monstertruck_modeling::{
    Face as MFace, Shell as MShell, Solid as MSolid, Surface as MSurface, builder, primitive,
};
use monstertruck_topology::shell::ShellCondition;

use super::{DraftError, draft};

const VOLUME_EPS: f64 = 1.0e-4;

/// Axis-aligned unit cube with planar faces and line edges.
fn unit_cube() -> MSolid {
    let bounding_box = BoundingBox::from_iter([Point3::origin(), Point3::new(1.0, 1.0, 1.0)]);
    primitive::cuboid(bounding_box)
}

/// Finds the id of the cube face whose outward normal is `normal`.
fn face_id_with_normal(solid: &MSolid, normal: Vector3) -> monstertruck_topology::FaceId<MSurface> {
    solid.boundaries()[0]
        .iter()
        .find(|face| face.oriented_surface().normal(0.5, 0.5).near(&normal))
        .map(|face| face.id())
        .unwrap()
}

fn tessellated_volume(solid: &MSolid) -> f64 { solid.triangulation(0.01).to_polygon().volume() }

/// The bottom face's plane (`z = 0`), used as the neutral plane in every
/// test.
fn bottom_plane() -> Plane { Plane::xy() }

/// Axis-aligned unit cube built via an extrude-chain (`vertex` -> `extrude`
/// along x -> `extrude` along y -> `extrude` along z), mirroring the
/// construction `openshape`'s `build_cuboid` uses.
///
/// Unlike [`unit_cube`] (`primitive::cuboid`), whose six faces are each
/// constructed directly from an explicitly hand-ordered `Plane::new(p, q, r)`
/// so that every face's `orientation()` is `true`, this construction routes
/// two of the side faces (the ones swept, by the third extrude, from the
/// "closing" edges the *second* extrude introduced to join the swept edge to
/// its translated copy) through `ExtrudeConnector`, which builds each side
/// face's plane from the swept edge's own direction and the sweep vector in
/// a fixed order regardless of which way is "outward". For two of the four
/// side faces that raw order happens to be inward, so the topology layer
/// marks those faces `orientation() == false` and corrects the *reported*
/// (`oriented_surface()`) normal by swapping the plane's `u`/`v` axes at
/// query time -- this is the shape that used to reproduce the
/// `drafted_shell` plane-orientation bug fixed in this module.
fn unit_cube_extrude_chain() -> MSolid {
    let v = builder::vertex(Point3::origin());
    let e = builder::extrude(&v, Vector3::unit_x());
    let f = builder::extrude(&e, Vector3::unit_y());
    builder::extrude(&f, Vector3::unit_z())
}

#[test]
fn draft_all_sides_gives_frustum_volume() {
    let cube = unit_cube();
    let sides = [
        Vector3::unit_x(),
        -Vector3::unit_x(),
        Vector3::unit_y(),
        -Vector3::unit_y(),
    ]
    .map(|normal| face_id_with_normal(&cube, normal));

    let drafted = draft(&cube, &sides, bottom_plane(), Vector3::unit_z(), Deg(5.0)).unwrap();

    assert_eq!(drafted.boundaries().len(), 1);
    assert_eq!(drafted.boundaries()[0].len(), 6);
    assert_eq!(
        drafted.boundaries()[0].shell_condition(),
        ShellCondition::Closed
    );

    // Each side moves inward by tan(5°) at the top; the base is untouched.
    let top_side = 1.0 - 2.0 * 5.0_f64.to_radians().tan();
    let (area_bottom, area_top) = (1.0, top_side * top_side);
    let expected = (area_bottom + area_top + (area_bottom * area_top).sqrt()) / 3.0;
    let volume = tessellated_volume(&drafted);
    assert!(
        (volume - expected).abs() < VOLUME_EPS,
        "volume {volume}, expected {expected}"
    );
}

/// Regression test for a fork bug: two of `unit_cube_extrude_chain`'s four
/// side faces have `orientation() == false` (see that function's doc
/// comment). `drafted_shell` used to rebuild every drafted face's surface
/// via `plane_through`, which always returns a plane whose *own* normal is
/// already the correct outward direction (it is built directly from
/// `face_plane`'s already-oriented-and-rotated normal), but then paired that
/// always-correct plane with the face's *unflipped* `absolute_boundaries()`
/// wire and re-applied the original face's orientation bit on top. For a
/// `false`-orientation face that re-applied flip inverted the (already
/// correct) new plane's effective normal instead of correcting the wire,
/// silently drafting that face's surface inward while its edges/vertices
/// were still moved correctly outward -- the shell stayed a valid closed
/// manifold (so nothing errored), but its enclosed volume came out wrong
/// (`0.9125...` instead of the frustum's `0.8352...`) because the
/// inward-facing face's triangles contributed a wrong-signed term to the
/// mesh volume integral.
#[test]
fn draft_all_sides_gives_frustum_volume_extrude_chain_cuboid() {
    let cube = unit_cube_extrude_chain();
    assert!(
        cube.boundaries()[0].iter().any(|face| !face.orientation()),
        "test cube must contain at least one orientation()==false face to exercise the bug"
    );
    let sides = [
        Vector3::unit_x(),
        -Vector3::unit_x(),
        Vector3::unit_y(),
        -Vector3::unit_y(),
    ]
    .map(|normal| face_id_with_normal(&cube, normal));

    let drafted = draft(&cube, &sides, bottom_plane(), Vector3::unit_z(), Deg(5.0)).unwrap();

    assert_eq!(drafted.boundaries().len(), 1);
    assert_eq!(drafted.boundaries()[0].len(), 6);
    assert_eq!(
        drafted.boundaries()[0].shell_condition(),
        ShellCondition::Closed
    );

    // Same frustum as `draft_all_sides_gives_frustum_volume`: each side
    // moves inward by tan(5°) at the top, the base is untouched.
    let top_side = 1.0 - 2.0 * 5.0_f64.to_radians().tan();
    let (area_bottom, area_top) = (1.0, top_side * top_side);
    let expected = (area_bottom + area_top + (area_bottom * area_top).sqrt()) / 3.0;
    let volume = tessellated_volume(&drafted);
    assert!(
        (volume - expected).abs() < VOLUME_EPS,
        "volume {volume}, expected {expected}"
    );
}

#[test]
fn draft_single_face_gives_wedge_volume() {
    let cube = unit_cube();
    let front = face_id_with_normal(&cube, -Vector3::unit_y());

    let drafted = draft(&cube, &[front], bottom_plane(), Vector3::unit_z(), Deg(5.0)).unwrap();

    assert_eq!(drafted.boundaries().len(), 1);
    assert_eq!(drafted.boundaries()[0].len(), 6);
    assert_eq!(
        drafted.boundaries()[0].shell_condition(),
        ShellCondition::Closed
    );

    // The front wall tilts inward by tan(5°) at the top, everything else is
    // fixed: a right-trapezoidal prism of unit depth, cross-section area
    // 1 - tan(5°)/2.
    let expected = 1.0 - 0.5 * 5.0_f64.to_radians().tan();
    let volume = tessellated_volume(&drafted);
    assert!(
        (volume - expected).abs() < VOLUME_EPS,
        "volume {volume}, expected {expected}"
    );
}

#[test]
fn draft_rejects_face_parallel_to_neutral_plane() {
    let cube = unit_cube();
    let top = face_id_with_normal(&cube, Vector3::unit_z());
    assert!(matches!(
        draft(&cube, &[top], bottom_plane(), Vector3::unit_z(), Deg(5.0)),
        Err(DraftError::ParallelToNeutralPlane)
    ));
}

#[test]
fn draft_rejects_excessive_angle() {
    let cube = unit_cube();
    let sides = [
        Vector3::unit_x(),
        -Vector3::unit_x(),
        Vector3::unit_y(),
        -Vector3::unit_y(),
    ]
    .map(|normal| face_id_with_normal(&cube, normal));
    // tan(60°) ≈ 1.73: opposing sides swing well past the cube's midline,
    // inverting the top face's edges.
    assert!(matches!(
        draft(&cube, &sides, bottom_plane(), Vector3::unit_z(), Deg(60.0)),
        Err(DraftError::ExcessiveAngle { .. })
    ));
}

#[test]
fn draft_rejects_curved_surface() {
    // A geometrically flat but NURBS-typed surface swapped in for the
    // cube's top face: planarity is decided by type, so this must be
    // rejected, not silently accepted.
    let cube = unit_cube();
    let shell = cube.boundaries()[0].clone();
    let top_index = shell
        .iter()
        .position(|face| {
            face.oriented_surface()
                .normal(0.5, 0.5)
                .near(&Vector3::unit_z())
        })
        .unwrap();
    let surface = BsplineSurface::new(
        (KnotVector::bezier_knot(1), KnotVector::bezier_knot(1)),
        vec![
            vec![Point3::new(0.0, 0.0, 1.0), Point3::new(0.0, 1.0, 1.0)],
            vec![Point3::new(1.0, 0.0, 1.0), Point3::new(1.0, 1.0, 1.0)],
        ],
    );
    let curved = MFace::new(
        shell[top_index].boundaries(),
        MSurface::BsplineSurface(surface),
    );
    let mixed: MShell = shell
        .iter()
        .enumerate()
        .map(|(index, face)| match index == top_index {
            true => curved.clone(),
            false => face.clone(),
        })
        .collect();
    let solid = MSolid::try_new(vec![mixed]).unwrap();
    let ids: Vec<_> = solid.boundaries()[0].iter().map(|face| face.id()).collect();
    assert!(matches!(
        draft(&solid, &ids, bottom_plane(), Vector3::unit_z(), Deg(5.0)),
        Err(DraftError::UnsupportedGeometry { .. })
    ));
}

#[test]
fn draft_rejects_invalid_angle() {
    let cube = unit_cube();
    let front = face_id_with_normal(&cube, -Vector3::unit_y());
    assert!(matches!(
        draft(&cube, &[front], bottom_plane(), Vector3::unit_z(), Deg(0.0)),
        Err(DraftError::InvalidAngle { .. })
    ));
    assert!(matches!(
        draft(
            &cube,
            &[front],
            bottom_plane(),
            Vector3::unit_z(),
            Deg(-5.0)
        ),
        Err(DraftError::InvalidAngle { .. })
    ));
    assert!(matches!(
        draft(
            &cube,
            &[front],
            bottom_plane(),
            Vector3::unit_z(),
            Deg(90.0)
        ),
        Err(DraftError::InvalidAngle { .. })
    ));
}

#[test]
fn draft_rejects_unknown_face() {
    let cube = unit_cube();
    let other = unit_cube();
    let foreign = face_id_with_normal(&other, -Vector3::unit_y());
    assert!(matches!(
        draft(
            &cube,
            &[foreign],
            bottom_plane(),
            Vector3::unit_z(),
            Deg(5.0)
        ),
        Err(DraftError::FaceNotFound)
    ));
}

#[test]
fn draft_rejects_no_faces_selected() {
    let cube = unit_cube();
    assert!(matches!(
        draft(&cube, &[], bottom_plane(), Vector3::unit_z(), Deg(5.0)),
        Err(DraftError::NoFacesSelected)
    ));
}
