//! Integration tests for [`path_sweep::try_path_sweep`].

use monstertruck_meshing::prelude::*;
use monstertruck_modeling::errors::Error;
use monstertruck_modeling::*;
use std::ops::Bound;

const VOLUME_EPS: f64 = 1.0e-3;

/// Builds a CCW square wire of side `2 * half` centered at `center`, spanning
/// the plane through `center` with in-plane axes `u` and `v`.
fn square_wire(center: Point3, u: Vector3, v: Vector3, half: f64) -> Wire {
    let v0 = builder::vertex(center - half * u - half * v);
    let v1 = builder::vertex(center + half * u - half * v);
    let v2 = builder::vertex(center + half * u + half * v);
    let v3 = builder::vertex(center - half * u + half * v);
    vec![
        builder::line(&v0, &v1),
        builder::line(&v1, &v2),
        builder::line(&v2, &v3),
        builder::line(&v3, &v0),
    ]
    .into()
}

/// Tessellated volume of a solid; also asserts the boundary shell is closed.
fn closed_volume(solid: &Solid) -> f64 {
    assert!(
        solid
            .boundaries()
            .iter()
            .all(|shell| shell.shell_condition() == shell::ShellCondition::Closed),
        "sweep result is not a closed shell"
    );
    solid.triangulation(0.01).to_polygon().volume()
}

/// A quarter-turn circular arc of radius `radius`, centered at the origin, in
/// the `y = 0` plane: `t in [0, pi/2] -> (radius cos t, 0, radius sin t)`.
#[derive(Clone, Copy)]
struct QuarterArc {
    radius: f64,
}

impl ParametricCurve for QuarterArc {
    type Point = Point3;
    type Vector = Vector3;
    fn evaluate(&self, t: f64) -> Point3 {
        Point3::new(self.radius * t.cos(), 0.0, self.radius * t.sin())
    }
    fn derivative(&self, t: f64) -> Vector3 {
        Vector3::new(-self.radius * t.sin(), 0.0, self.radius * t.cos())
    }
    fn derivative_2(&self, t: f64) -> Vector3 {
        Vector3::new(-self.radius * t.cos(), 0.0, -self.radius * t.sin())
    }
    fn derivative_n(&self, n: usize, t: f64) -> Vector3 {
        match n % 4 {
            0 => Vector3::new(self.radius * t.cos(), 0.0, self.radius * t.sin()),
            1 => self.derivative(t),
            2 => self.derivative_2(t),
            _ => Vector3::new(self.radius * t.sin(), 0.0, -self.radius * t.cos()),
        }
    }
    fn parameter_range(&self) -> ParameterRange {
        (
            Bound::Included(0.0),
            Bound::Included(std::f64::consts::FRAC_PI_2),
        )
    }
}
impl BoundedCurve for QuarterArc {}

/// A helix-like S-curve: `t in [0, t_end] -> (r cos t, pitch * t, r sin t)`.
#[derive(Clone, Copy)]
struct Helix {
    radius: f64,
    pitch: f64,
    t_end: f64,
}

impl ParametricCurve for Helix {
    type Point = Point3;
    type Vector = Vector3;
    fn evaluate(&self, t: f64) -> Point3 {
        Point3::new(self.radius * t.cos(), self.pitch * t, self.radius * t.sin())
    }
    fn derivative(&self, t: f64) -> Vector3 {
        Vector3::new(-self.radius * t.sin(), self.pitch, self.radius * t.cos())
    }
    fn derivative_2(&self, t: f64) -> Vector3 {
        Vector3::new(-self.radius * t.cos(), 0.0, -self.radius * t.sin())
    }
    fn derivative_n(&self, n: usize, t: f64) -> Vector3 {
        match n {
            0 => self.evaluate(t).to_vec(),
            1 => self.derivative(t),
            2 => self.derivative_2(t),
            _ => Vector3::new(self.radius * t.sin(), 0.0, -self.radius * t.cos()),
        }
    }
    fn parameter_range(&self) -> ParameterRange {
        (Bound::Included(0.0), Bound::Included(self.t_end))
    }
}
impl BoundedCurve for Helix {}

/// A curve that reports a period, standing in for any closed/periodic path.
#[derive(Clone, Copy)]
struct PeriodicStub;

impl ParametricCurve for PeriodicStub {
    type Point = Point3;
    type Vector = Vector3;
    fn evaluate(&self, t: f64) -> Point3 { Point3::new(t.cos(), t.sin(), 0.0) }
    fn derivative(&self, t: f64) -> Vector3 { Vector3::new(-t.sin(), t.cos(), 0.0) }
    fn derivative_2(&self, t: f64) -> Vector3 { Vector3::new(-t.cos(), -t.sin(), 0.0) }
    fn derivative_n(&self, _n: usize, t: f64) -> Vector3 { self.derivative(t) }
    fn parameter_range(&self) -> ParameterRange {
        (
            Bound::Included(0.0),
            Bound::Included(2.0 * std::f64::consts::PI),
        )
    }
    fn period(&self) -> Option<f64> { Some(2.0 * std::f64::consts::PI) }
}
impl BoundedCurve for PeriodicStub {}

// -- (a) straight-line path ≡ extrude --

#[test]
fn straight_line_matches_extrude_volume() {
    let square = square_wire(Point3::origin(), Vector3::unit_x(), Vector3::unit_y(), 0.5);
    let path = Line(Point3::new(0.0, 0.0, 0.0), Point3::new(0.0, 0.0, 1.0));
    let solid: Solid = path_sweep::try_path_sweep(&square, &path, 5).unwrap();
    let volume = closed_volume(&solid);
    assert!(
        (volume - 1.0).abs() < VOLUME_EPS,
        "expected unit volume, got {volume}"
    );
}

#[test]
fn straight_line_matches_extrude_solid() {
    let square = square_wire(Point3::origin(), Vector3::unit_x(), Vector3::unit_y(), 0.5);
    let path = Line(Point3::new(0.0, 0.0, 0.0), Point3::new(0.0, 0.0, 1.0));
    let swept: Solid = path_sweep::try_path_sweep(&square, &path, 4).unwrap();

    let disk: Face = builder::try_attach_plane(std::slice::from_ref(&square)).unwrap();
    let extruded: Solid = builder::extrude(&disk, Vector3::new(0.0, 0.0, 1.0));

    let swept_volume = closed_volume(&swept);
    let extruded_volume = closed_volume(&extruded);
    assert!(
        (swept_volume - extruded_volume).abs() < VOLUME_EPS,
        "swept volume {swept_volume} != extruded volume {extruded_volume}"
    );
}

// -- (b) quarter-circle arc path, square profile ≈ Pappus' theorem --

#[test]
fn quarter_arc_matches_pappus_theorem() {
    let radius = 3.0;
    let half = 0.4; // well clear of the axis (the world y axis) at x = radius.
    let square = square_wire(
        Point3::new(radius, 0.0, 0.0),
        Vector3::unit_x(),
        Vector3::unit_y(),
        half,
    );
    let path = QuarterArc { radius };
    let solid: Solid = path_sweep::try_path_sweep(&square, &path, 12).unwrap();
    let volume = closed_volume(&solid);

    // Pappus: V = area * (distance traveled by the centroid).
    // The centroid sits at distance `radius` from the sweep axis (the world
    // y axis), and travels a quarter of that circle's circumference.
    let area = (2.0 * half) * (2.0 * half);
    let centroid_path_length = radius * std::f64::consts::FRAC_PI_2;
    let expected = area * centroid_path_length;
    assert!(
        (volume - expected).abs() < 1.0e-2,
        "expected Pappus volume {expected}, got {volume}"
    );
}

// -- (c) helix-ish S-curve path: closed, positive-volume shell --

#[test]
fn helix_path_produces_closed_positive_volume_solid() {
    let helix = Helix {
        radius: 2.0,
        pitch: 0.3,
        t_end: 4.0,
    };
    let tangent0 = helix.derivative(0.0).normalize();
    // Build an in-plane basis transverse to the start tangent so the profile
    // starts out perpendicular to the path, per `try_path_sweep`'s documented
    // convention.
    let hint = if tangent0.x.abs() < 0.9 {
        Vector3::unit_x()
    } else {
        Vector3::unit_y()
    };
    let u = tangent0.cross(hint).normalize();
    let v = tangent0.cross(u).normalize();
    let square = square_wire(helix.evaluate(0.0), u, v, 0.2);

    let solid: Solid = path_sweep::try_path_sweep(&square, &helix, 24).unwrap();
    let volume = closed_volume(&solid);
    assert!(volume > 0.0, "expected positive volume, got {volume}");
}

// -- (d) degenerate inputs --

#[test]
fn zero_divisions_is_an_error() {
    let square = square_wire(Point3::origin(), Vector3::unit_x(), Vector3::unit_y(), 0.5);
    let path = Line(Point3::new(0.0, 0.0, 0.0), Point3::new(0.0, 0.0, 1.0));
    let result: Result<Solid> = path_sweep::try_path_sweep(&square, &path, 0);
    assert_eq!(result.unwrap_err(), Error::TooFewDivisions);
}

#[test]
fn one_division_succeeds() {
    let square = square_wire(Point3::origin(), Vector3::unit_x(), Vector3::unit_y(), 0.5);
    let path = Line(Point3::new(0.0, 0.0, 0.0), Point3::new(0.0, 0.0, 1.0));
    let solid: Solid = path_sweep::try_path_sweep(&square, &path, 1).unwrap();
    let volume = closed_volume(&solid);
    assert!(
        (volume - 1.0).abs() < VOLUME_EPS,
        "expected unit volume, got {volume}"
    );
}

#[test]
fn open_profile_wire_is_an_error() {
    let v0 = builder::vertex(Point3::new(-0.5, -0.5, 0.0));
    let v1 = builder::vertex(Point3::new(0.5, -0.5, 0.0));
    let v2 = builder::vertex(Point3::new(0.5, 0.5, 0.0));
    let open: Wire = vec![builder::line(&v0, &v1), builder::line(&v1, &v2)].into();
    let path = Line(Point3::new(0.0, 0.0, 0.0), Point3::new(0.0, 0.0, 1.0));
    let result: Result<Solid> = path_sweep::try_path_sweep(&open, &path, 4);
    assert_eq!(result.unwrap_err(), Error::OpenWire);
}

#[test]
fn periodic_path_is_rejected() {
    let square = square_wire(
        Point3::new(1.0, 0.0, 0.0),
        Vector3::unit_z(),
        Vector3::unit_y(),
        0.1,
    );
    let result: Result<Solid> = path_sweep::try_path_sweep(&square, &PeriodicStub, 8);
    assert_eq!(result.unwrap_err(), Error::ClosedPathNotSupported);
}
