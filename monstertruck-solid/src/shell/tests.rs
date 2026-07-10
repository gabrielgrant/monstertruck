use monstertruck_geometry::prelude::*;
use monstertruck_meshing::prelude::*;
use monstertruck_modeling::{
    Curve as MCurve, Edge as MEdge, Face as MFace, Shell as MShell, Solid as MSolid,
    Surface as MSurface, Vertex as MVertex, Wire as MWire, primitive,
};
use monstertruck_topology::shell::ShellCondition;

use super::{ShellError, hollow, thicken};

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

/// A rectangular planar face spanning `min..max` in the coordinates of `plane`.
fn rect_face(min: Point2, max: Point2, plane: Plane) -> MFace {
    let wire: MWire = primitive::rect(BoundingBox::from_iter([min, max]), plane);
    MFace::new(vec![wire], MSurface::Plane(plane))
}

fn tessellated_volume(solid: &MSolid) -> f64 { solid.triangulation(0.01).to_polygon().volume() }

#[test]
fn hollow_unit_cube_open_top() {
    let cube = unit_cube();
    let top = face_id_with_normal(&cube, Vector3::unit_z());

    let hollowed = hollow(&cube, 0.1, &[top]).unwrap();

    // 5 kept outer faces + 1 ring + 5 inner cavity faces.
    assert_eq!(hollowed.boundaries().len(), 1);
    assert_eq!(hollowed.boundaries()[0].len(), 11);
    assert_eq!(
        hollowed.boundaries()[0].shell_condition(),
        ShellCondition::Closed
    );

    // The cavity spans [0.1, 0.9] x [0.1, 0.9] x [0.1, 1.0]: it reaches the
    // open face's plane, so its height is 0.9, not 0.8.
    let expected = 1.0 - 0.8 * 0.8 * 0.9;
    let volume = tessellated_volume(&hollowed);
    assert!(
        (volume - expected).abs() < VOLUME_EPS,
        "volume {volume}, expected {expected}"
    );
}

#[test]
fn hollow_unit_cube_closed() {
    let cube = unit_cube();

    let hollowed = hollow(&cube, 0.1, &[]).unwrap();

    // Outer boundary and inner cavity boundary.
    assert_eq!(hollowed.boundaries().len(), 2);
    for boundary in hollowed.boundaries() {
        assert_eq!(boundary.shell_condition(), ShellCondition::Closed);
    }

    let expected = 1.0 - 0.8 * 0.8 * 0.8;
    let volume = tessellated_volume(&hollowed);
    assert!(
        (volume - expected).abs() < VOLUME_EPS,
        "volume {volume}, expected {expected}"
    );
}

#[test]
fn thicken_flat_sheet() {
    let sheet: MShell = [rect_face(
        Point2::new(0.0, 0.0),
        Point2::new(1.0, 2.0),
        Plane::xy(),
    )]
    .into();

    let slab = thicken(&sheet, 0.25).unwrap();

    // Original face + offset face + 4 side walls.
    assert_eq!(slab.boundaries().len(), 1);
    assert_eq!(slab.boundaries()[0].len(), 6);
    assert_eq!(
        slab.boundaries()[0].shell_condition(),
        ShellCondition::Closed
    );

    let expected = 1.0 * 2.0 * 0.25;
    let volume = tessellated_volume(&slab);
    assert!(
        (volume - expected).abs() < VOLUME_EPS,
        "volume {volume}, expected {expected}"
    );
}

#[test]
fn thicken_open_box() {
    // A unit cube with its top face removed: an open shell whose boundary
    // wire is the top rim. Thickening it is equivalent to hollowing the
    // cube with the top open.
    let cube = unit_cube();
    let top = face_id_with_normal(&cube, Vector3::unit_z());
    let open_box: MShell = cube.boundaries()[0]
        .iter()
        .filter(|face| face.id() != top)
        .cloned()
        .collect();
    assert_eq!(open_box.shell_condition(), ShellCondition::Oriented);

    let thickened = thicken(&open_box, 0.1).unwrap();

    // 5 outer faces + 5 inner faces + 4 rim walls.
    assert_eq!(thickened.boundaries().len(), 1);
    assert_eq!(thickened.boundaries()[0].len(), 14);

    let expected = 1.0 - 0.8 * 0.8 * 0.9;
    let volume = tessellated_volume(&thickened);
    assert!(
        (volume - expected).abs() < VOLUME_EPS,
        "volume {volume}, expected {expected}"
    );
}

#[test]
fn thicken_closed_shell_gives_hollow_solid() {
    let cube = unit_cube();

    let hollowed = thicken(&cube.boundaries()[0], 0.1).unwrap();

    assert_eq!(hollowed.boundaries().len(), 2);
    let expected = 1.0 - 0.8 * 0.8 * 0.8;
    let volume = tessellated_volume(&hollowed);
    assert!(
        (volume - expected).abs() < VOLUME_EPS,
        "volume {volume}, expected {expected}"
    );
}

#[test]
fn thicken_rejects_nonpositive_thickness() {
    let sheet: MShell = [rect_face(
        Point2::new(0.0, 0.0),
        Point2::new(1.0, 1.0),
        Plane::xy(),
    )]
    .into();
    assert!(matches!(
        thicken(&sheet, 0.0),
        Err(ShellError::NonPositiveThickness(_))
    ));
    assert!(matches!(
        thicken(&sheet, -0.5),
        Err(ShellError::NonPositiveThickness(_))
    ));
}

#[test]
fn hollow_rejects_excessive_thickness() {
    // 0.6 > half the cube width: the cavity edges would invert.
    let cube = unit_cube();
    assert!(matches!(
        hollow(&cube, 0.6, &[]),
        Err(ShellError::ThicknessTooLarge { .. })
    ));
    // Exactly half the width: the cavity collapses to a point.
    assert!(matches!(
        hollow(&cube, 0.5, &[]),
        Err(ShellError::ThicknessTooLarge { .. })
    ));
}

#[test]
fn hollow_rejects_unknown_face() {
    let cube = unit_cube();
    let other = unit_cube();
    let foreign = face_id_with_normal(&other, Vector3::unit_z());
    assert!(matches!(
        hollow(&cube, 0.1, &[foreign]),
        Err(ShellError::FaceNotFound)
    ));
}

#[test]
fn hollow_rejects_adjacent_open_faces() {
    let cube = unit_cube();
    let top = face_id_with_normal(&cube, Vector3::unit_z());
    let front = face_id_with_normal(&cube, -Vector3::unit_y());
    assert!(matches!(
        hollow(&cube, 0.1, &[top, front]),
        Err(ShellError::AdjacentOpenFaces)
    ));
}

#[test]
fn hollow_rejects_already_hollow_solid() {
    let cube = unit_cube();
    let hollowed = hollow(&cube, 0.1, &[]).unwrap();
    assert!(matches!(
        hollow(&hollowed, 0.01, &[]),
        Err(ShellError::UnsupportedSolidBoundaries(2))
    ));
}

#[test]
fn thicken_rejects_curved_surface() {
    // A geometrically flat but NURBS-typed surface: planarity is decided by
    // type, so this must be rejected, not silently accepted.
    let p = [
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(1.0, 1.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
    ];
    let v = MVertex::from_points(p);
    let line = |i: usize, j: usize| MEdge::new(&v[i], &v[j], MCurve::Line(Line(p[i], p[j])));
    let wire: MWire = [line(0, 1), line(1, 2), line(2, 3), line(3, 0)].into();
    let surface = BsplineSurface::new(
        (KnotVector::bezier_knot(1), KnotVector::bezier_knot(1)),
        vec![vec![p[0], p[3]], vec![p[1], p[2]]],
    );
    let face = MFace::new(vec![wire], MSurface::BsplineSurface(surface));
    let shell: MShell = [face].into();
    assert!(matches!(
        thicken(&shell, 0.1),
        Err(ShellError::UnsupportedGeometry { .. })
    ));
}

#[test]
fn thicken_rejects_disconnected_shell() {
    let far_plane = Plane::new(
        Point3::new(5.0, 0.0, 0.0),
        Point3::new(6.0, 0.0, 0.0),
        Point3::new(5.0, 1.0, 0.0),
    );
    let shell: MShell = [
        rect_face(Point2::new(0.0, 0.0), Point2::new(1.0, 1.0), Plane::xy()),
        rect_face(Point2::new(0.0, 0.0), Point2::new(1.0, 1.0), far_plane),
    ]
    .into();
    assert!(matches!(
        thicken(&shell, 0.1),
        Err(ShellError::NotOrientable { .. })
    ));
}
