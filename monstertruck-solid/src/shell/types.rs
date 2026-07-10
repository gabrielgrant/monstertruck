use monstertruck_geometry::prelude::*;

monstertruck_topology::prelude!(Point3, Line<Point3>, Plane, pub(super));

/// Curve types that can participate in shell/thicken operations: anything
/// convertible losslessly to and from a straight [`Line`].
///
/// Automatically implemented for any type satisfying the bounds.
pub trait ThickenableCurve:
    Clone + ParametricCurve<Point = Point3> + BoundedCurve + TryInto<Line<Point3>> + From<Line<Point3>>
{
}

impl<T> ThickenableCurve for T where T: Clone
        + ParametricCurve<Point = Point3>
        + BoundedCurve
        + TryInto<Line<Point3>>
        + From<Line<Point3>>
{
}

/// Surface types that can participate in shell/thicken operations: anything
/// convertible losslessly to and from a [`Plane`].
///
/// Automatically implemented for any type satisfying the bounds.
pub trait ThickenableSurface:
    Clone + ParametricSurface<Point = Point3> + TryInto<Plane> + From<Plane> {
}

impl<T> ThickenableSurface for T where T: Clone + ParametricSurface<Point = Point3> + TryInto<Plane> + From<Plane> {}

/// Converts an external shell into the internal planar ([`Line`]/[`Plane`])
/// representation. Returns `None` if any surface is not a plane or any edge
/// curve is not a straight line.
pub(super) fn convert_shell_in<C: ThickenableCurve, S: ThickenableSurface>(
    shell: &monstertruck_topology::Shell<Point3, C, S>,
) -> Option<Shell> {
    shell.try_mapped(
        |p| Some(*p),
        |c| c.clone().try_into().ok(),
        |s| s.clone().try_into().ok(),
    )
}

/// Converts an internal planar solid back to the caller's curve/surface
/// representation.
pub(super) fn convert_solid_out<C: ThickenableCurve, S: ThickenableSurface>(
    solid: &Solid,
) -> monstertruck_topology::Solid<Point3, C, S> {
    solid.mapped(
        |p| *p,
        |c: &Line<Point3>| C::from(*c),
        |s: &Plane| S::from(*s),
    )
}
