//! Sweeps a closed profile wire along an arbitrary 3D path curve, using a
//! rotation-minimizing frame to avoid the twisting artifacts of a naive
//! Frenet frame.
//!
//! See [`try_path_sweep`] for the entry point.

use crate::{Result, errors::Error};
use monstertruck_geometry::prelude::*;

type Wire<C> = monstertruck_topology::Wire<Point3, C>;
type Shell<C, S> = monstertruck_topology::Shell<Point3, C, S>;
type Solid<C, S> = monstertruck_topology::Solid<Point3, C, S>;

/// A single rotation-minimizing sweep frame: a point on the path together
/// with an orthonormal right-handed triad `(r, tangent.cross(r), tangent)`.
#[derive(Clone, Copy, Debug)]
struct Frame {
    origin: Point3,
    tangent: Vector3,
    r: Vector3,
}

/// Picks a unit vector orthogonal to `n`, deterministically, by crossing `n`
/// with whichever world axis it is least parallel to.
///
/// This only seeds the very first frame; the rotation-minimizing property
/// makes the sweep's shape independent of this choice (it only fixes the
/// "phase" of the frame around the tangent, i.e. how the profile's local x
/// axis is initially oriented).
fn arbitrary_perpendicular(n: Vector3) -> Vector3 {
    let a = n.map(f64::abs);
    if a.x > a.z || a.y > a.z {
        Vector3::new(-n.y, n.x, 0.0).normalize()
    } else {
        Vector3::new(-n.z, 0.0, n.x).normalize()
    }
}

/// Advances a rotation-minimizing frame from `(x0, t0, r0)` to the frame at
/// `(x1, t1)`, using the double-reflection method of Wang, Jüttler, Zheng &
/// Liu, *"Computation of Rotation Minimizing Frames"*, ACM TOG 27(1), 2008.
///
/// Unlike a Frenet frame, this stays well-defined on straight segments and
/// at inflection points, since it never depends on the curve's second
/// derivative.
fn double_reflection(x0: Point3, t0: Vector3, r0: Vector3, x1: Point3, t1: Vector3) -> Vector3 {
    let v1 = x1 - x0;
    let c1 = v1.dot(v1);
    let (r_l, t_l) = if c1.so_small() {
        // `x1` coincides with `x0` (a degenerate, zero-length step): the
        // positional reflection is undefined, so skip it and reflect only
        // across the tangent change below.
        (r0, t0)
    } else {
        (
            r0 - (2.0 / c1) * v1.dot(r0) * v1,
            t0 - (2.0 / c1) * v1.dot(t0) * v1,
        )
    };
    let v2 = t1 - t_l;
    let c2 = v2.dot(v2);
    if c2.so_small() {
        r_l.normalize()
    } else {
        (r_l - (2.0 / c2) * v2.dot(r_l) * v2).normalize()
    }
}

/// Samples `path` at parameters `ts` and propagates a rotation-minimizing
/// frame across the samples.
///
/// # Errors
/// [`Error::DegenerateSweepTangent`] if the path's tangent vanishes at any
/// sample parameter.
fn rotation_minimizing_frames<C: ParametricCurve3D>(path: &C, ts: &[f64]) -> Result<Vec<Frame>> {
    let mut frames: Vec<Frame> = Vec::with_capacity(ts.len());
    for &t in ts {
        let origin = path.evaluate(t);
        let tangent_raw = path.derivative(t);
        if tangent_raw.magnitude2().so_small() {
            return Err(Error::DegenerateSweepTangent);
        }
        let tangent = tangent_raw.normalize();
        let r = match frames.last() {
            None => arbitrary_perpendicular(tangent),
            Some(prev) if tangent.near(&prev.tangent) => {
                // Exact handling of the straight (constant-tangent) case:
                // skip the reflection algebra entirely so a dead-straight
                // path accumulates zero drift from floating-point error.
                prev.r
            }
            Some(prev) => double_reflection(prev.origin, prev.tangent, prev.r, origin, tangent),
        };
        frames.push(Frame { origin, tangent, r });
    }
    Ok(frames)
}

/// Returns the rigid transform (rotation + translation) carrying the frame
/// `base` onto `frame`.
fn frame_transform(base: &Frame, frame: &Frame) -> Matrix4 {
    let base_rot = Matrix3::from_cols(base.r, base.tangent.cross(base.r), base.tangent);
    let rot = Matrix3::from_cols(frame.r, frame.tangent.cross(frame.r), frame.tangent);
    // `base_rot` is orthonormal, so its inverse is its transpose.
    let m = rot * base_rot.transpose();
    let translation = frame.origin.to_vec() - m * base.origin.to_vec();
    Matrix4::from_translation(translation) * Matrix4::from(m)
}

/// Sweeps a closed `profile` wire along `path`, producing a solid tube.
///
/// The path is sampled at `divisions + 1` evenly spaced parameters and the
/// profile is rigidly transported to each sample using a rotation-minimizing
/// frame (the double-reflection method; see [`double_reflection`]), rather
/// than a naive Frenet frame, so the sweep does not twist unpredictably on
/// straight segments or through inflection points. The resulting cross
/// sections are skinned together with [`try_skin_wires`](crate::builder::try_skin_wires)
/// and the two ends are capped with [`try_attach_plane`](crate::builder::try_attach_plane).
///
/// # Convention
///
/// `profile` is used exactly as given for the first cross section -- it is
/// *not* re-centered or reprojected. Every later cross section is obtained
/// by applying the rigid motion that carries the path's start frame to the
/// frame at that sample, pivoting about `path.evaluate(t0)`. In other words:
///
/// ```text
/// cross_section(t) = path(t) + R(t) * (profile - path(t0))
/// ```
///
/// where `R(t)` is the rotation-minimizing rotation taking the tangent at
/// `t0` to the tangent at `t`. For a sensible tube, `profile` should lie in
/// (or near) the plane through `path.evaluate(t0)` perpendicular to
/// `path.derivative(t0)` -- e.g. a unit square in the `z = 0` plane swept
/// along the `z` axis reproduces [`extrude`](crate::builder::extrude).
///
/// # Errors
///
/// - [`Error::OpenWire`] if `profile` is not closed.
/// - [`Error::TooFewDivisions`] if `divisions == 0`.
/// - [`Error::DegenerateSweepTangent`] if `path`'s tangent vanishes at some
///   sample.
/// - [`Error::ClosedPathNotSupported`] if `path` reports a `period()` (a
///   periodic/closed path); closing the tube into a torus-like solid instead
///   of capping it is not yet implemented.
/// - Propagates errors from `try_skin_wires`, `try_attach_plane`, and
///   `Solid::try_new`.
///
/// # Examples
///
/// ```
/// use monstertruck_modeling::*;
///
/// // A unit square swept one unit along the z axis is a unit cube.
/// let v0 = builder::vertex(Point3::new(-0.5, -0.5, 0.0));
/// let v1 = builder::vertex(Point3::new(0.5, -0.5, 0.0));
/// let v2 = builder::vertex(Point3::new(0.5, 0.5, 0.0));
/// let v3 = builder::vertex(Point3::new(-0.5, 0.5, 0.0));
/// let square: Wire = vec![
///     builder::line(&v0, &v1),
///     builder::line(&v1, &v2),
///     builder::line(&v2, &v3),
///     builder::line(&v3, &v0),
/// ]
/// .into();
/// let path = Line(Point3::new(0.0, 0.0, 0.0), Point3::new(0.0, 0.0, 1.0));
/// let solid: Solid = path_sweep::try_path_sweep(&square, &path, 3).unwrap();
/// assert_eq!(solid.boundaries()[0].shell_condition(), shell::ShellCondition::Closed);
/// ```
pub fn try_path_sweep<Cp, Cw, S>(
    profile: &Wire<Cp>,
    path: &Cw,
    divisions: usize,
) -> Result<Solid<Cp, S>>
where
    Cp: ParametricCurve3D + BoundedCurve + Invertible + Transformed<Matrix4>,
    Cw: ParametricCurve3D + BoundedCurve,
    Line<Point3>: ToSameGeometry<Cp>,
    HomotopySurface<Cp, Cp>: ToSameGeometry<S>,
    Plane: IncludeCurve<Cp> + ToSameGeometry<S>,
{
    if profile.is_empty() {
        return Err(Error::FromTopology(
            monstertruck_topology::errors::Error::EmptyWire,
        ));
    }
    if !profile.is_closed() {
        return Err(Error::OpenWire);
    }
    if divisions < 1 {
        return Err(Error::TooFewDivisions);
    }
    if path.period().is_some() {
        return Err(Error::ClosedPathNotSupported);
    }

    let (t0, t1) = path.range_tuple();
    let ts: Vec<f64> = (0..=divisions)
        .map(|i| t0 + (t1 - t0) * (i as f64 / divisions as f64))
        .collect();
    let frames = rotation_minimizing_frames(path, &ts)?;

    let base = frames[0];
    let wires: Vec<Wire<Cp>> = frames
        .iter()
        .map(|frame| crate::builder::transformed(profile, frame_transform(&base, frame)))
        .collect();

    let mut shell: Shell<Cp, S> = crate::builder::try_skin_wires(&wires)?;
    let start_cap = crate::builder::try_attach_plane(vec![wires[0].clone()])?.inverse();
    let end_cap = crate::builder::try_attach_plane(vec![wires[wires.len() - 1].clone()])?;
    shell.push(start_cap);
    shell.push(end_cap);

    Ok(Solid::try_new(vec![shell])?)
}
