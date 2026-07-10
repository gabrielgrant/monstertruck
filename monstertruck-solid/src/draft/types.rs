use monstertruck_geometry::prelude::*;

monstertruck_topology::prelude!(Point3, Line<Point3>, Plane, pub(super));

/// Curve types that can participate in draft operations: anything
/// convertible losslessly to and from a straight [`Line`].
///
/// Automatically implemented for any type satisfying the bounds. Identical
/// in shape to the shell/thicken module's `ThickenableCurve`; kept as a
/// separate trait here because the two modules do not share a dependency
/// edge on this branch (see the crate-level notes in the draft module).
pub trait DraftableCurve:
    Clone + ParametricCurve<Point = Point3> + BoundedCurve + TryInto<Line<Point3>> + From<Line<Point3>>
{
}

impl<T> DraftableCurve for T where T: Clone
        + ParametricCurve<Point = Point3>
        + BoundedCurve
        + TryInto<Line<Point3>>
        + From<Line<Point3>>
{
}

/// Surface types that can participate in draft operations: anything
/// convertible losslessly to and from a [`Plane`].
///
/// Automatically implemented for any type satisfying the bounds.
pub trait DraftableSurface:
    Clone + ParametricSurface<Point = Point3> + TryInto<Plane> + From<Plane> {
}

impl<T> DraftableSurface for T where T: Clone + ParametricSurface<Point = Point3> + TryInto<Plane> + From<Plane> {}

/// Converts an external shell into the internal planar ([`Line`]/[`Plane`])
/// representation. Returns `None` if any surface is not a plane or any edge
/// curve is not a straight line.
pub(super) fn convert_shell_in<C: DraftableCurve, S: DraftableSurface>(
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
pub(super) fn convert_solid_out<C: DraftableCurve, S: DraftableSurface>(
    solid: &Solid,
) -> monstertruck_topology::Solid<Point3, C, S> {
    solid.mapped(
        |p| *p,
        |c: &Line<Point3>| C::from(*c),
        |s: &Plane| S::from(*s),
    )
}
