use super::*;
use derive_more::From;
use monstertruck_geometry::prelude::{
    BoundaryCurve2D, SupportsExactPatchDomains, TryIntoBsplineSurface,
    TryIntoHomogeneousBsplineCurve, TryIntoHomogeneousBsplineSurface,
};
#[doc(hidden)]
pub use monstertruck_geometry::prelude::{algo, inv_or_zero};
pub use monstertruck_geometry::{decorators::*, nurbs::*, specifieds::*, t_spline::*};
pub use monstertruck_mesh::PolylineCurve;
use monstertruck_topology::compress::{CompressedTrimmedShell, CompressedTrimmedSolid};
use monstertruck_topology::trimmed::{TrimmedFace, TrimmedShell, TrimmedSolid};
use monstertruck_traits::SnapCurveEndpoints;
#[cfg(not(target_arch = "wasm32"))]
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

type ModelSurfaceCurve = SurfaceCurve<
    Box<Curve>,
    Box<Surface>,
    Box<Surface>,
    ParameterCurve<Curve2D, Box<Surface>>,
    ParameterCurve<Curve2D, Box<Surface>>,
>;

/// 3-dimensional curve
#[derive(
    Clone,
    Debug,
    Serialize,
    Deserialize,
    From,
    ParametricCurve,
    BoundedCurve,
    Cut,
    Invertible,
    SearchNearestParameterD1,
    SearchParameterD1,
)]
pub enum Curve {
    /// line
    Line(Line<Point3>),
    /// 3-dimensional B-spline curve
    BsplineCurve(BsplineCurve<Point3>),
    /// 3-dimensional NURBS curve
    NurbsCurve(NurbsCurve<Vector4>),
    /// 3-dimensional curve carried by a 2-dimensional parameter curve on a surface
    #[allow(clippy::enum_variant_names)]
    ParameterCurve(ParameterCurve<Curve2D, Box<Surface>>),
    /// intersection curve
    IntersectionCurve(ModelSurfaceCurve),
}

/// 2-dimensional curve used as a parameter-space trim.
#[derive(
    Clone,
    Debug,
    Serialize,
    Deserialize,
    From,
    ParametricCurve,
    BoundedCurve,
    ParameterDivision1D,
    Cut,
    Invertible,
    SearchNearestParameterD1,
    SearchParameterD1,
)]
pub enum Conic2D {
    /// ellipse
    Ellipse(Processor<TrimmedCurve<UnitCircle<Point2>>, Matrix3>),
    /// hyperbola
    Hyperbola(Processor<TrimmedCurve<UnitHyperbola<Point2>>, Matrix3>),
    /// parabola
    Parabola(Processor<TrimmedCurve<UnitParabola<Point2>>, Matrix3>),
}

/// 2-dimensional curve used as a parameter-space trim.
#[derive(
    Clone,
    Debug,
    Serialize,
    Deserialize,
    From,
    ParametricCurve,
    BoundedCurve,
    ParameterDivision1D,
    Cut,
    Invertible,
    SearchNearestParameterD1,
    SearchParameterD1,
)]
pub enum Curve2D {
    /// line
    Line(Line<Point2>),
    /// polyline
    Polyline(PolylineCurve<Point2>),
    /// conic
    Conic(Conic2D),
    /// 2-dimensional B-spline curve
    BsplineCurve(BsplineCurve<Point2>),
    /// 2-dimensional NURBS curve
    NurbsCurve(NurbsCurve<Vector3>),
}

macro_rules! derive_curve_method {
    ($curve: expr, $method: expr, $($ver: ident),*) => {
        match $curve {
            Curve::Line(got) => $method(got, $($ver), *),
            Curve::BsplineCurve(got) => $method(got, $($ver), *),
            Curve::NurbsCurve(got) => $method(got, $($ver), *),
            Curve::ParameterCurve(got) => $method(got, $($ver), *),
            Curve::IntersectionCurve(got) => $method(got, $($ver), *),
        }
    };
}

macro_rules! derive_curve_self_method {
    ($curve: expr, $method: expr, $($ver: ident),*) => {
        match $curve {
            Curve::Line(got) => Curve::Line($method(got, $($ver), *)),
            Curve::BsplineCurve(got) => Curve::BsplineCurve($method(got, $($ver), *)),
            Curve::NurbsCurve(got) => Curve::NurbsCurve($method(got, $($ver), *)),
            Curve::ParameterCurve(got) => Curve::ParameterCurve($method(got, $($ver), *)),
            Curve::IntersectionCurve(got) => Curve::IntersectionCurve($method(got, $($ver), *)),
        }
    };
}

fn sample_curve_to_nurbs(curve: &(impl ParametricCurve3D + BoundedCurve)) -> NurbsCurve<Vector4> {
    let (t0, t1) = curve.range_tuple();
    let samples = 16usize;
    let points: Vec<Point3> = (0..=samples)
        .map(|i| t0 + (t1 - t0) * (i as f64) / (samples as f64))
        .map(|t| curve.evaluate(t))
        .collect();
    let knots: Vec<f64> = (0..=samples).map(|i| i as f64 / samples as f64).collect();
    let knot_vec = KnotVector::from(
        std::iter::once(0.0)
            .chain(knots.iter().copied())
            .chain(std::iter::once(1.0))
            .collect::<Vec<_>>(),
    );
    NurbsCurve::from(BsplineCurve::new(knot_vec, points))
}

fn linear_bspline_division(
    curve: &BsplineCurve<Point3>,
    range: (f64, f64),
) -> Option<(Vec<f64>, Vec<Point3>)> {
    let curve_range = curve.range_tuple();
    (curve.degree() == 1 && curve_range.0.near(&range.0) && curve_range.1.near(&range.1)).then(
        || {
            (
                (1..=curve.control_points().len())
                    .map(|index| curve.knot(index))
                    .collect(),
                curve.control_points().clone(),
            )
        },
    )
}

fn sampled_parameter_boundary(
    curve: &(
         impl ParametricCurve3D<Point = Point3> + BoundedCurve + ParameterDivision1D<Point = Point3>
     ),
    surface: &Surface,
    tolerance: f64,
) -> Option<Vec<Point2>> {
    let points = curve.parameter_division(curve.range_tuple(), tolerance).1;
    let project = |point: Point3, hint: Option<(f64, f64)>| {
        surface
            .search_nearest_parameter(point, hint, 100)
            .or_else(|| surface.search_parameter(point, hint, 100))
            .or_else(|| surface.search_nearest_parameter(point, None, 100))
            .or_else(|| surface.search_parameter(point, None, 100))
            .map(|(u, v)| Point2::new(u, v))
    };
    points
        .iter()
        .copied()
        .scan(None, |hint, point| {
            let uv = project(point, *hint);
            *hint = uv.map(|uv| (uv.x, uv.y));
            Some(uv)
        })
        .collect::<Option<Vec<_>>>()
        .or_else(|| {
            points
                .into_iter()
                .map(|point| project(point, None))
                .collect()
        })
}

fn curve2d_from_sampled_boundary(points: Vec<Point2>) -> Option<Curve2D> {
    if points.len() < 2 {
        None
    } else {
        let front = points.first().copied()?;
        let back = points.last().copied()?;
        let line = Line(front, back);
        let is_linear = points.iter().copied().all(|point| {
            line.search_nearest_parameter(point, None, 1)
                .is_some_and(|t| line.evaluate(t).near(&point))
        });
        if is_linear {
            Some(Curve2D::Line(line))
        } else {
            let denom = (points.len() - 1) as f64;
            let knot_vec = KnotVector::from(
                std::iter::once(0.0)
                    .chain((0..points.len()).map(|index| index as f64 / denom))
                    .chain(std::iter::once(1.0))
                    .collect::<Vec<_>>(),
            );
            Some(Curve2D::BsplineCurve(BsplineCurve::new(knot_vec, points)))
        }
    }
}

fn same_surface(lhs: &Surface, rhs: &Surface) -> bool {
    if std::mem::discriminant(lhs) != std::mem::discriminant(rhs) {
        false
    } else if let (Some((lu0, lu1)), Some((lv0, lv1)), Some((ru0, ru1)), Some((rv0, rv1))) = (
        lhs.try_range_tuple().0,
        lhs.try_range_tuple().1,
        rhs.try_range_tuple().0,
        rhs.try_range_tuple().1,
    ) {
        [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0), (1.0, 1.0), (0.5, 0.5)]
            .into_iter()
            .all(|(s, t)| {
                let lp = lhs.evaluate(lu0 + (lu1 - lu0) * s, lv0 + (lv1 - lv0) * t);
                let rp = rhs.evaluate(ru0 + (ru1 - ru0) * s, rv0 + (rv1 - rv0) * t);
                lp.near(&rp)
            })
    } else {
        false
    }
}

fn exact_line_boundary(
    line: &Line<Point3>,
    surface: &Surface,
) -> Option<ParameterCurve<Curve2D, Box<Surface>>> {
    match surface {
        Surface::Plane(plane) => line.exact_parameter_boundary_2d(plane).map(|boundary| {
            let (curve, plane) = boundary.decompose();
            ParameterCurve::new(Curve2D::Line(curve), Box::new(Surface::Plane(plane)))
        }),
        Surface::RevolutionSurface(surface) => {
            line.exact_parameter_boundary_2d(surface).map(|boundary| {
                let (curve, surface) = boundary.decompose();
                ParameterCurve::new(
                    Curve2D::Line(curve),
                    Box::new(Surface::RevolutionSurface(surface)),
                )
            })
        }
        _ => None,
    }
}

fn exact_bspline_boundary(
    curve: &BsplineCurve<Point3>,
    surface: &Surface,
) -> Option<ParameterCurve<Curve2D, Box<Surface>>> {
    match surface {
        Surface::BsplineSurface(surface) => {
            curve.exact_parameter_boundary_2d(surface).map(|boundary| {
                let (curve, surface) = boundary.decompose();
                ParameterCurve::new(
                    Curve2D::Line(curve),
                    Box::new(Surface::BsplineSurface(surface)),
                )
            })
        }
        Surface::RevolutionSurface(surface) => {
            curve.exact_parameter_boundary_2d(surface).map(|boundary| {
                let (curve, surface) = boundary.decompose();
                ParameterCurve::new(
                    Curve2D::Line(curve),
                    Box::new(Surface::RevolutionSurface(surface)),
                )
            })
        }
        _ => None,
    }
}

fn exact_nurbs_boundary(
    curve: &NurbsCurve<Vector4>,
    surface: &Surface,
) -> Option<ParameterCurve<Curve2D, Box<Surface>>> {
    match surface {
        Surface::NurbsSurface(surface) => {
            curve.exact_parameter_boundary_2d(surface).map(|boundary| {
                let (curve, surface) = boundary.decompose();
                ParameterCurve::new(
                    Curve2D::Line(curve),
                    Box::new(Surface::NurbsSurface(surface)),
                )
            })
        }
        Surface::RevolutionSurface(surface) => {
            curve.exact_parameter_boundary_2d(surface).map(|boundary| {
                let (curve, surface) = boundary.decompose();
                ParameterCurve::new(
                    Curve2D::Line(curve),
                    Box::new(Surface::RevolutionSurface(surface)),
                )
            })
        }
        _ => None,
    }
}

impl Transformed<Matrix4> for Curve {
    fn transform_by(&mut self, trans: Matrix4) {
        derive_curve_method!(self, Transformed::transform_by, trans);
    }
    fn transformed(&self, trans: Matrix4) -> Self {
        derive_curve_self_method!(self, Transformed::transformed, trans)
    }
}

impl ParameterDivision1D for Curve {
    type Point = Point3;
    fn parameter_division(&self, range: (f64, f64), tol: f64) -> (Vec<f64>, Vec<Self::Point>) {
        let debug_profile = std::env::var("MT_PROFILE_CURVE_DIVISION").is_ok();
        let started = std::time::Instant::now();
        let result = match self {
            Curve::Line(curve) => curve.parameter_division(range, tol),
            Curve::BsplineCurve(curve) => linear_bspline_division(curve, range)
                .unwrap_or_else(|| curve.parameter_division(range, tol)),
            Curve::NurbsCurve(curve) => curve.parameter_division(range, tol),
            Curve::ParameterCurve(curve) => curve.parameter_division(range, tol),
            Curve::IntersectionCurve(curve) => curve.leader().parameter_division(range, tol),
        };
        if debug_profile {
            let kind = match self {
                Curve::Line(_) => "Line",
                Curve::BsplineCurve(_) => "BsplineCurve",
                Curve::NurbsCurve(_) => "NurbsCurve",
                Curve::ParameterCurve(_) => "ParameterCurve",
                Curve::IntersectionCurve(curve) => match curve.leader().as_ref() {
                    Curve::Line(_) => "IntersectionCurve(Line)",
                    Curve::BsplineCurve(_) => "IntersectionCurve(BsplineCurve)",
                    Curve::NurbsCurve(_) => "IntersectionCurve(NurbsCurve)",
                    Curve::ParameterCurve(_) => "IntersectionCurve(ParameterCurve)",
                    Curve::IntersectionCurve(_) => "IntersectionCurve(IntersectionCurve)",
                },
            };
            eprintln!(
                "trace bool model_curve_division kind={} points={} tol={} elapsed_ms={:.3}",
                kind,
                result.1.len(),
                tol,
                started.elapsed().as_secs_f64() * 1000.0,
            );
        }
        result
    }
}

impl From<IntersectionCurve<BsplineCurve<Point3>, Surface, Surface>> for Curve {
    fn from(c: IntersectionCurve<BsplineCurve<Point3>, Surface, Surface>) -> Curve {
        let (surface0, surface1, leader) = c.destruct();
        Curve::IntersectionCurve(SurfaceCurve::with_boundaries(
            Box::new(surface0),
            Box::new(surface1),
            Box::new(leader.into()),
            None,
            None,
        ))
    }
}

fn boundary_curve_2d_to_model_curve(
    curve: ParameterCurve<BoundaryCurve2D, Surface>,
) -> ParameterCurve<Curve2D, Box<Surface>> {
    let surface = Box::new(curve.surface().clone());
    match curve.curve() {
        BoundaryCurve2D::Line(line) => ParameterCurve::new(Curve2D::Line(*line), surface),
        BoundaryCurve2D::BsplineCurve(bspline) => {
            ParameterCurve::new(Curve2D::BsplineCurve(bspline.clone()), surface)
        }
        BoundaryCurve2D::NurbsCurve(nurbs) => {
            ParameterCurve::new(Curve2D::NurbsCurve(nurbs.clone()), surface)
        }
    }
}

impl
    From<
        SurfaceCurve<
            BsplineCurve<Point3>,
            Surface,
            Surface,
            ParameterCurve<BoundaryCurve2D, Surface>,
            ParameterCurve<BoundaryCurve2D, Surface>,
        >,
    > for Curve
{
    fn from(
        c: SurfaceCurve<
            BsplineCurve<Point3>,
            Surface,
            Surface,
            ParameterCurve<BoundaryCurve2D, Surface>,
            ParameterCurve<BoundaryCurve2D, Surface>,
        >,
    ) -> Curve {
        let (surface0, surface1, leader, boundary0, boundary1) = c.destruct_with_boundaries();
        Curve::IntersectionCurve(SurfaceCurve::with_boundaries(
            Box::new(surface0),
            Box::new(surface1),
            Box::new(leader.into()),
            boundary0.map(boundary_curve_2d_to_model_curve),
            boundary1.map(boundary_curve_2d_to_model_curve),
        ))
    }
}

impl
    From<
        SurfaceCurve<
            NurbsCurve<Vector4>,
            Surface,
            Surface,
            ParameterCurve<BoundaryCurve2D, Surface>,
            ParameterCurve<BoundaryCurve2D, Surface>,
        >,
    > for Curve
{
    fn from(
        c: SurfaceCurve<
            NurbsCurve<Vector4>,
            Surface,
            Surface,
            ParameterCurve<BoundaryCurve2D, Surface>,
            ParameterCurve<BoundaryCurve2D, Surface>,
        >,
    ) -> Curve {
        let (surface0, surface1, leader, boundary0, boundary1) = c.destruct_with_boundaries();
        Curve::IntersectionCurve(SurfaceCurve::with_boundaries(
            Box::new(surface0),
            Box::new(surface1),
            Box::new(leader.into()),
            boundary0.map(boundary_curve_2d_to_model_curve),
            boundary1.map(boundary_curve_2d_to_model_curve),
        ))
    }
}

impl ToSameGeometry<Curve> for Line<Point3> {
    #[inline]
    fn to_same_geometry(&self) -> Curve { Curve::from(*self) }
}

impl ToSameGeometry<Curve> for Processor<TrimmedCurve<UnitCircle<Point3>>, Matrix4> {
    #[inline]
    fn to_same_geometry(&self) -> Curve { Curve::NurbsCurve(self.to_same_geometry()) }
}

// Lossless extraction of a straight line, used by planar-only algorithms
// such as `monstertruck_solid::shell` and `monstertruck_solid::draft`.
impl TryFrom<Curve> for Line<Point3> {
    type Error = ();
    fn try_from(curve: Curve) -> std::result::Result<Self, ()> {
        match curve {
            Curve::Line(line) => Ok(line),
            _ => Err(()),
        }
    }
}

impl ToSameGeometry<Curve> for BsplineCurve<Point3> {
    #[inline]
    fn to_same_geometry(&self) -> Curve { Curve::from(self.clone()) }
}

impl Curve {
    /// Into non-ratinalized 4-dimensional B-spline curve
    pub fn lift_up(&self) -> BsplineCurve<Vector4> {
        match self {
            Curve::Line(curve) => Curve::BsplineCurve((*curve).into()).lift_up(),
            Curve::BsplineCurve(curve) => BsplineCurve::new(
                curve.knot_vector().clone(),
                curve
                    .control_points()
                    .iter()
                    .map(|pt| pt.to_vec().extend(1.0))
                    .collect(),
            ),
            Curve::NurbsCurve(curve) => curve.non_rationalized().clone(),
            Curve::ParameterCurve(curve) => sample_curve_to_nurbs(curve).non_rationalized().clone(),
            Curve::IntersectionCurve(curve) => curve.leader().lift_up(),
        }
    }
}

fn curve2d_endpoint_hints(curve: &Curve2D) -> Option<(Point2, Point2)> {
    match curve {
        Curve2D::Line(curve) => Some((curve.0, curve.1)),
        Curve2D::Polyline(curve) => Some((*curve.first()?, *curve.last()?)),
        Curve2D::BsplineCurve(curve) => Some((
            *curve.control_points().first()?,
            *curve.control_points().last()?,
        )),
        Curve2D::NurbsCurve(curve) => Some((
            curve.control_points().first()?.to_point(),
            curve.control_points().last()?.to_point(),
        )),
        Curve2D::Conic(_) => None,
    }
}

fn set_curve2d_endpoints(curve: &mut Curve2D, front: Point2, back: Point2) {
    match curve {
        Curve2D::Line(curve) => {
            curve.0 = front;
            curve.1 = back;
        }
        Curve2D::Polyline(curve) => {
            if let Some(point) = curve.first_mut() {
                *point = front;
            }
            if let Some(point) = curve.last_mut() {
                *point = back;
            }
        }
        Curve2D::BsplineCurve(curve) => {
            if !curve.control_points().is_empty() {
                *curve.control_point_mut(0) = front;
            }
            if curve.control_points().len() > 1 {
                let last = curve.control_points().len() - 1;
                *curve.control_point_mut(last) = back;
            }
        }
        Curve2D::NurbsCurve(curve) => {
            if !curve.control_points().is_empty() {
                let point = curve.control_point_mut(0);
                let weight = point.weight();
                *point = front.to_vec().extend(weight);
            }
            if curve.control_points().len() > 1 {
                let last = curve.control_points().len() - 1;
                let point = curve.control_point_mut(last);
                let weight = point.weight();
                *point = back.to_vec().extend(weight);
            }
        }
        Curve2D::Conic(_) => {}
    }
}

fn snap_parameter_curve_endpoints(
    curve: &mut ParameterCurve<Curve2D, Box<Surface>>,
    front: Point3,
    back: Point3,
) {
    let hints = curve2d_endpoint_hints(curve.curve());
    let front_uv = hints
        .and_then(|(front_hint, _)| {
            curve
                .surface()
                .search_nearest_parameter(front, Some((front_hint.x, front_hint.y)), 100)
        })
        .or_else(|| curve.surface().search_nearest_parameter(front, None, 100))
        .map(|(u, v)| Point2::new(u, v));
    let back_uv = hints
        .and_then(|(_, back_hint)| {
            curve
                .surface()
                .search_nearest_parameter(back, Some((back_hint.x, back_hint.y)), 100)
        })
        .or_else(|| curve.surface().search_nearest_parameter(back, None, 100))
        .map(|(u, v)| Point2::new(u, v));
    if let (Some(front_uv), Some(back_uv)) = (front_uv, back_uv) {
        let (mut curve2d, surface) = curve.clone().decompose();
        set_curve2d_endpoints(&mut curve2d, front_uv, back_uv);
        *curve = ParameterCurve::new(curve2d, surface);
    }
}

impl SnapCurveEndpoints for Curve {
    fn snap_endpoints(&mut self, front: Point3, back: Point3) {
        match self {
            Curve::IntersectionCurve(curve) => {
                curve.leader_mut().snap_endpoints(front, back);
                if let Some(boundary) = curve.boundary0_mut() {
                    snap_parameter_curve_endpoints(boundary, front, back);
                }
                if let Some(boundary) = curve.boundary1_mut() {
                    snap_parameter_curve_endpoints(boundary, front, back);
                }
            }
            Curve::ParameterCurve(curve) => snap_parameter_curve_endpoints(curve, front, back),
            _ => {}
        }
    }
}

impl TryIntoHomogeneousBsplineCurve for Curve {
    fn try_into_homogeneous_bspline_curve(&self) -> Option<BsplineCurve<Vector4>> {
        match self {
            Curve::Line(curve) => curve.try_into_homogeneous_bspline_curve(),
            Curve::BsplineCurve(curve) => curve.try_into_homogeneous_bspline_curve(),
            Curve::NurbsCurve(curve) => curve.try_into_homogeneous_bspline_curve(),
            Curve::ParameterCurve(_) => None,
            Curve::IntersectionCurve(curve) => curve.leader().try_into_homogeneous_bspline_curve(),
        }
    }
}

impl TryFrom<ParameterCurve<Line<Point2>, Surface>> for Curve {
    type Error = ();
    fn try_from(curve: ParameterCurve<Line<Point2>, Surface>) -> std::result::Result<Self, ()> {
        let (line, surface) = curve.decompose();
        Ok(Curve::ParameterCurve(ParameterCurve::new(
            Curve2D::Line(line),
            Box::new(surface),
        )))
    }
}

impl ParameterBoundary2D<Surface> for Curve {
    fn parameter_boundary_2d(&self, surface: &Surface, tolerance: f64) -> Option<Vec<Point2>> {
        match self {
            Curve::Line(curve) => {
                sampled_parameter_boundary(curve, surface, tolerance).or_else(|| {
                    exact_line_boundary(curve, surface).map(|boundary| {
                        boundary
                            .curve()
                            .parameter_division(boundary.curve().range_tuple(), tolerance)
                            .1
                    })
                })
            }
            Curve::BsplineCurve(curve) => sampled_parameter_boundary(curve, surface, tolerance)
                .or_else(|| {
                    exact_bspline_boundary(curve, surface).map(|boundary| {
                        boundary
                            .curve()
                            .parameter_division(boundary.curve().range_tuple(), tolerance)
                            .1
                    })
                }),
            Curve::NurbsCurve(curve) => sampled_parameter_boundary(curve, surface, tolerance)
                .or_else(|| {
                    exact_nurbs_boundary(curve, surface).map(|boundary| {
                        boundary
                            .curve()
                            .parameter_division(boundary.curve().range_tuple(), tolerance)
                            .1
                    })
                }),
            Curve::ParameterCurve(curve) => {
                same_surface(curve.surface().as_ref(), surface).then(|| {
                    curve
                        .curve()
                        .parameter_division(curve.curve().range_tuple(), tolerance)
                        .1
                })
            }
            Curve::IntersectionCurve(curve) => {
                if same_surface(curve.surface0().as_ref(), surface) {
                    curve
                        .boundary0()
                        .map(|boundary| {
                            boundary
                                .curve()
                                .parameter_division(boundary.curve().range_tuple(), tolerance)
                                .1
                        })
                        .or_else(|| curve.leader().parameter_boundary_2d(surface, tolerance))
                } else if same_surface(curve.surface1().as_ref(), surface) {
                    curve
                        .boundary1()
                        .map(|boundary| {
                            boundary
                                .curve()
                                .parameter_division(boundary.curve().range_tuple(), tolerance)
                                .1
                        })
                        .or_else(|| curve.leader().parameter_boundary_2d(surface, tolerance))
                } else {
                    sampled_parameter_boundary(curve.leader().as_ref(), surface, tolerance)
                }
            }
        }
    }
}

impl ExactParameterBoundary2D<Surface> for Curve {
    type BoundaryCurve = ParameterCurve<Curve2D, Box<Surface>>;

    fn exact_parameter_boundary_2d(&self, surface: &Surface) -> Option<Self::BoundaryCurve> {
        match self {
            Curve::Line(curve) => exact_line_boundary(curve, surface),
            Curve::BsplineCurve(curve) => exact_bspline_boundary(curve, surface),
            Curve::NurbsCurve(curve) => exact_nurbs_boundary(curve, surface),
            Curve::ParameterCurve(curve) if same_surface(curve.surface().as_ref(), surface) => {
                Some(curve.clone())
            }
            Curve::IntersectionCurve(curve) if same_surface(curve.surface0().as_ref(), surface) => {
                curve
                    .boundary0()
                    .cloned()
                    .or_else(|| curve.leader().exact_parameter_boundary_2d(surface))
            }
            Curve::IntersectionCurve(curve) if same_surface(curve.surface1().as_ref(), surface) => {
                curve
                    .boundary1()
                    .cloned()
                    .or_else(|| curve.leader().exact_parameter_boundary_2d(surface))
            }
            _ => None,
        }
    }
}

impl BoundaryCurveFromSamples<Surface> for ParameterCurve<Curve2D, Box<Surface>> {
    fn boundary_curve_from_samples(surface: &Surface, points: Vec<Point2>) -> Option<Self> {
        curve2d_from_sampled_boundary(points)
            .map(|curve| ParameterCurve::new(curve, Box::new(surface.clone())))
    }
}

impl Curve {
    /// Converts this curve into a face-local parameter curve on `surface`.
    ///
    /// Exact trim data is preserved when available. Otherwise this falls back
    /// to a sampled polyline trim in the surface domain.
    pub fn to_parameter_curve_on(
        &self,
        surface: &Surface,
        tolerance: f64,
    ) -> Option<ParameterCurve<Curve2D, Box<Surface>>> {
        let debug_profile = std::env::var("MT_PROFILE_PARAMETER_CURVE_ON").is_ok();
        let started = std::time::Instant::now();
        let exact = self.exact_parameter_boundary_2d(surface);
        let exact_hit = exact.is_some();
        let result = exact.or_else(|| {
            self.parameter_boundary_2d(surface, tolerance)
                .filter(|points| points.len() >= 2)
                .and_then(curve2d_from_sampled_boundary)
                .map(|curve| ParameterCurve::new(curve, Box::new(surface.clone())))
        });
        if debug_profile {
            let kind = match self {
                Curve::Line(_) => "Line",
                Curve::BsplineCurve(_) => "BsplineCurve",
                Curve::NurbsCurve(_) => "NurbsCurve",
                Curve::ParameterCurve(_) => "ParameterCurve",
                Curve::IntersectionCurve(curve) => match curve.leader().as_ref() {
                    Curve::Line(_) => "IntersectionCurve(Line)",
                    Curve::BsplineCurve(_) => "IntersectionCurve(BsplineCurve)",
                    Curve::NurbsCurve(_) => "IntersectionCurve(NurbsCurve)",
                    Curve::ParameterCurve(_) => "IntersectionCurve(ParameterCurve)",
                    Curve::IntersectionCurve(_) => "IntersectionCurve(IntersectionCurve)",
                },
            };
            eprintln!(
                "trace bool parameter_curve_on kind={} exact={} output={} elapsed_ms={:.3}",
                kind,
                exact_hit,
                result.is_some(),
                started.elapsed().as_secs_f64() * 1000.0,
            );
        }
        result
    }
}

/// Extension trait for creating runtime trimmed topology with face-local parameter curves.
pub trait ToTrimmedParameterCurves {
    /// The trimmed topology output.
    type Output;

    /// Creates runtime trimmed topology with face-local parameter curves.
    fn to_trimmed_with_parameter_curves(&self, tolerance: f64) -> Self::Output;
}

impl ToTrimmedParameterCurves for Shell {
    type Output = TrimmedShell<Point3, Curve, Surface, ParameterCurve<Curve2D, Box<Surface>>>;

    #[cfg(not(target_arch = "wasm32"))]
    fn to_trimmed_with_parameter_curves(&self, tolerance: f64) -> Self::Output {
        self.iter()
            .cloned()
            .collect::<Vec<_>>()
            .into_par_iter()
            .map(|face| {
                let surface = face.surface();
                let trims = face
                    .absolute_boundaries()
                    .iter()
                    .map(|wire| {
                        wire.iter()
                            .map(|edge| edge.curve().to_parameter_curve_on(&surface, tolerance))
                            .collect()
                    })
                    .collect();
                TrimmedFace::new(face, trims)
            })
            .collect::<Vec<_>>()
            .into_iter()
            .collect()
    }

    #[cfg(target_arch = "wasm32")]
    fn to_trimmed_with_parameter_curves(&self, tolerance: f64) -> Self::Output {
        self.to_trimmed_with_face_trims(|edge, surface| {
            edge.curve().to_parameter_curve_on(surface, tolerance)
        })
    }
}

impl ToTrimmedParameterCurves for Solid {
    type Output = TrimmedSolid<Point3, Curve, Surface, ParameterCurve<Curve2D, Box<Surface>>>;

    #[cfg(not(target_arch = "wasm32"))]
    fn to_trimmed_with_parameter_curves(&self, tolerance: f64) -> Self::Output {
        TrimmedSolid::new(
            self.boundaries()
                .par_iter()
                .map(|shell| shell.to_trimmed_with_parameter_curves(tolerance))
                .collect(),
        )
    }

    #[cfg(target_arch = "wasm32")]
    fn to_trimmed_with_parameter_curves(&self, tolerance: f64) -> Self::Output {
        self.to_trimmed_with_face_trims(|edge, surface| {
            edge.curve().to_parameter_curve_on(surface, tolerance)
        })
    }
}

/// Extension trait for creating compressed trimmed topology with face-local parameter curves.
pub trait ToCompressedTrimmedParameterCurves {
    /// The compressed trimmed topology output.
    type Output;

    /// Creates compressed trimmed topology with face-local parameter curves.
    fn compress_with_parameter_curves(&self, tolerance: f64) -> Self::Output;
}

impl ToCompressedTrimmedParameterCurves for Shell {
    type Output =
        CompressedTrimmedShell<Point3, Curve, Surface, ParameterCurve<Curve2D, Box<Surface>>>;

    fn compress_with_parameter_curves(&self, tolerance: f64) -> Self::Output {
        let trimmed = self.to_trimmed_with_parameter_curves(tolerance);
        CompressedTrimmedShell::from(&trimmed)
    }
}

impl ToCompressedTrimmedParameterCurves for Solid {
    type Output =
        CompressedTrimmedSolid<Point3, Curve, Surface, ParameterCurve<Curve2D, Box<Surface>>>;

    fn compress_with_parameter_curves(&self, tolerance: f64) -> Self::Output {
        let trimmed = self.to_trimmed_with_parameter_curves(tolerance);
        CompressedTrimmedSolid::from(&trimmed)
    }
}

/// 3-dimensional surfaces
#[derive(
    Clone,
    Debug,
    Serialize,
    Deserialize,
    From,
    ParametricSurface,
    ParameterDivision2D,
    Invertible,
    SearchParameterD2,
)]
#[allow(clippy::large_enum_variant)]
pub enum Surface {
    /// Plane
    Plane(Plane),
    /// 3-dimensional B-spline surface
    BsplineSurface(BsplineSurface<Point3>),
    /// 3-dimensional NURBS Surface
    NurbsSurface(NurbsSurface<Vector4>),
    /// revoluted curve
    #[serde(alias = "RevolutedCurve")]
    RevolutionSurface(Processor<RevolutionSurface<Curve>, Matrix4>),
    /// T-spline surface
    TsplineSurface(Tmesh<Point3>),
}

macro_rules! derive_surface_method {
    ($surface: expr, $method: expr, $($ver: ident),*) => {
        match $surface {
            Self::Plane(got) => $method(got, $($ver), *),
            Self::BsplineSurface(got) => $method(got, $($ver), *),
            Self::NurbsSurface(got) => $method(got, $($ver), *),
            Self::RevolutionSurface(got) => $method(got, $($ver), *),
            Self::TsplineSurface(got) => $method(got, $($ver), *),
        }
    };
}

macro_rules! derive_surface_self_method {
    ($surface: expr, $method: expr, $($ver: ident),*) => {
        match $surface {
            Self::Plane(got) => Self::Plane($method(got, $($ver), *)),
            Self::BsplineSurface(got) => Self::BsplineSurface($method(got, $($ver), *)),
            Self::NurbsSurface(got) => Self::NurbsSurface($method(got, $($ver), *)),
            Self::RevolutionSurface(got) => Self::RevolutionSurface($method(got, $($ver), *)),
            Self::TsplineSurface(got) => Self::TsplineSurface($method(got, $($ver), *)),
        }
    };
}

impl ParametricSurface3D for Surface {
    #[inline(always)]
    fn normal(&self, u: f64, v: f64) -> Vector3 {
        derive_surface_method!(self, ParametricSurface3D::normal, u, v)
    }
}

impl Transformed<Matrix4> for Surface {
    fn transform_by(&mut self, trans: Matrix4) {
        derive_surface_method!(self, Transformed::transform_by, trans);
    }
    fn transformed(&self, trans: Matrix4) -> Self {
        derive_surface_self_method!(self, Transformed::transformed, trans)
    }
}

impl IncludeCurve<Curve> for Surface {
    #[inline(always)]
    fn include(&self, curve: &Curve) -> bool {
        if let Curve::ParameterCurve(curve) = curve {
            same_surface(curve.surface().as_ref(), self)
        } else {
            match self {
                Surface::BsplineSurface(surface) => match curve {
                    &Curve::Line(curve) => surface.include(&BsplineCurve::from(curve)),
                    Curve::BsplineCurve(curve) => surface.include(curve),
                    Curve::NurbsCurve(curve) => surface.include(curve),
                    Curve::IntersectionCurve(curve) => match curve.leader().as_ref() {
                        Curve::Line(curve) => surface.include(&BsplineCurve::from(*curve)),
                        Curve::BsplineCurve(curve) => surface.include(curve),
                        Curve::NurbsCurve(curve) => surface.include(curve),
                        Curve::ParameterCurve(_) => false,
                        Curve::IntersectionCurve(_) => false,
                    },
                    Curve::ParameterCurve(_) => unreachable!(),
                },
                Surface::NurbsSurface(surface) => match curve {
                    &Curve::Line(curve) => surface.include(&BsplineCurve::from(curve)),
                    Curve::BsplineCurve(curve) => surface.include(curve),
                    Curve::NurbsCurve(curve) => surface.include(curve),
                    Curve::IntersectionCurve(curve) => match curve.leader().as_ref() {
                        Curve::Line(curve) => surface.include(&BsplineCurve::from(*curve)),
                        Curve::BsplineCurve(curve) => surface.include(curve),
                        Curve::NurbsCurve(curve) => surface.include(curve),
                        Curve::ParameterCurve(_) => false,
                        Curve::IntersectionCurve(_) => false,
                    },
                    Curve::ParameterCurve(_) => unreachable!(),
                },
                Surface::Plane(surface) => match curve {
                    &Curve::Line(curve) => surface.include(&BsplineCurve::from(curve)),
                    Curve::BsplineCurve(curve) => surface.include(curve),
                    Curve::NurbsCurve(curve) => surface.include(curve),
                    Curve::IntersectionCurve(curve) => match curve.leader().as_ref() {
                        Curve::Line(curve) => surface.include(&BsplineCurve::from(*curve)),
                        Curve::BsplineCurve(curve) => surface.include(curve),
                        Curve::NurbsCurve(curve) => surface.include(curve),
                        Curve::ParameterCurve(_) => false,
                        Curve::IntersectionCurve(_) => false,
                    },
                    Curve::ParameterCurve(_) => unreachable!(),
                },
                Surface::TsplineSurface(surface) => {
                    curve.lift_up().control_points().iter().all(|v| {
                        let p = v.to_point();
                        surface.search_parameter(p, None, 1).is_some()
                    })
                }
                Surface::RevolutionSurface(surface) => match surface.entity_curve() {
                    &Curve::Line(entity_line) => {
                        let entity_bsp = BsplineCurve::from(entity_line);
                        let surface = RevolutionSurface::by_revolution(
                            &entity_bsp,
                            surface.origin(),
                            surface.axis(),
                        );
                        match curve {
                            &Curve::Line(curve) => surface.include(&BsplineCurve::from(curve)),
                            Curve::BsplineCurve(curve) => surface.include(curve),
                            Curve::NurbsCurve(curve) => surface.include(curve),
                            Curve::IntersectionCurve(curve) => match curve.leader().as_ref() {
                                Curve::Line(curve) => surface.include(&BsplineCurve::from(*curve)),
                                Curve::BsplineCurve(curve) => surface.include(curve),
                                Curve::NurbsCurve(curve) => surface.include(curve),
                                Curve::ParameterCurve(_) => false,
                                Curve::IntersectionCurve(_) => false,
                            },
                            Curve::ParameterCurve(_) => unreachable!(),
                        }
                    }
                    Curve::BsplineCurve(entity_curve) => {
                        let surface = RevolutionSurface::by_revolution(
                            entity_curve,
                            surface.origin(),
                            surface.axis(),
                        );
                        match curve {
                            &Curve::Line(curve) => surface.include(&BsplineCurve::from(curve)),
                            Curve::BsplineCurve(curve) => surface.include(curve),
                            Curve::NurbsCurve(curve) => surface.include(curve),
                            Curve::IntersectionCurve(curve) => match curve.leader().as_ref() {
                                Curve::Line(curve) => surface.include(&BsplineCurve::from(*curve)),
                                Curve::BsplineCurve(curve) => surface.include(curve),
                                Curve::NurbsCurve(curve) => surface.include(curve),
                                Curve::ParameterCurve(_) => false,
                                Curve::IntersectionCurve(_) => false,
                            },
                            Curve::ParameterCurve(_) => unreachable!(),
                        }
                    }
                    Curve::NurbsCurve(entity_curve) => {
                        let surface = RevolutionSurface::by_revolution(
                            entity_curve,
                            surface.origin(),
                            surface.axis(),
                        );
                        match curve {
                            &Curve::Line(curve) => surface.include(&BsplineCurve::from(curve)),
                            Curve::BsplineCurve(curve) => surface.include(curve),
                            Curve::NurbsCurve(curve) => surface.include(curve),
                            Curve::IntersectionCurve(curve) => match curve.leader().as_ref() {
                                Curve::Line(curve) => surface.include(&BsplineCurve::from(*curve)),
                                Curve::BsplineCurve(curve) => surface.include(curve),
                                Curve::NurbsCurve(curve) => surface.include(curve),
                                Curve::ParameterCurve(_) => false,
                                Curve::IntersectionCurve(_) => false,
                            },
                            Curve::ParameterCurve(_) => unreachable!(),
                        }
                    }
                    Curve::IntersectionCurve(entity_curve) => {
                        let leader = entity_curve.leader().as_ref();
                        match leader {
                            Curve::Line(entity_line) => {
                                let entity_bsp = BsplineCurve::from(*entity_line);
                                let surface = RevolutionSurface::by_revolution(
                                    &entity_bsp,
                                    surface.origin(),
                                    surface.axis(),
                                );
                                match curve {
                                    &Curve::Line(curve) => {
                                        surface.include(&BsplineCurve::from(curve))
                                    }
                                    Curve::BsplineCurve(curve) => surface.include(curve),
                                    Curve::NurbsCurve(curve) => surface.include(curve),
                                    Curve::IntersectionCurve(curve) => {
                                        match curve.leader().as_ref() {
                                            Curve::Line(curve) => {
                                                surface.include(&BsplineCurve::from(*curve))
                                            }
                                            Curve::BsplineCurve(curve) => surface.include(curve),
                                            Curve::NurbsCurve(curve) => surface.include(curve),
                                            Curve::ParameterCurve(_) => false,
                                            Curve::IntersectionCurve(_) => false,
                                        }
                                    }
                                    Curve::ParameterCurve(_) => unreachable!(),
                                }
                            }
                            Curve::BsplineCurve(entity_curve) => {
                                let surface = RevolutionSurface::by_revolution(
                                    entity_curve,
                                    surface.origin(),
                                    surface.axis(),
                                );
                                match curve {
                                    &Curve::Line(curve) => {
                                        surface.include(&BsplineCurve::from(curve))
                                    }
                                    Curve::BsplineCurve(curve) => surface.include(curve),
                                    Curve::NurbsCurve(curve) => surface.include(curve),
                                    Curve::IntersectionCurve(curve) => {
                                        match curve.leader().as_ref() {
                                            Curve::Line(curve) => {
                                                surface.include(&BsplineCurve::from(*curve))
                                            }
                                            Curve::BsplineCurve(curve) => surface.include(curve),
                                            Curve::NurbsCurve(curve) => surface.include(curve),
                                            Curve::ParameterCurve(_) => false,
                                            Curve::IntersectionCurve(_) => false,
                                        }
                                    }
                                    Curve::ParameterCurve(_) => unreachable!(),
                                }
                            }
                            Curve::NurbsCurve(entity_curve) => {
                                let surface = RevolutionSurface::by_revolution(
                                    entity_curve,
                                    surface.origin(),
                                    surface.axis(),
                                );
                                match curve {
                                    &Curve::Line(curve) => {
                                        surface.include(&BsplineCurve::from(curve))
                                    }
                                    Curve::BsplineCurve(curve) => surface.include(curve),
                                    Curve::NurbsCurve(curve) => surface.include(curve),
                                    Curve::IntersectionCurve(curve) => {
                                        match curve.leader().as_ref() {
                                            Curve::Line(curve) => {
                                                surface.include(&BsplineCurve::from(*curve))
                                            }
                                            Curve::BsplineCurve(curve) => surface.include(curve),
                                            Curve::NurbsCurve(curve) => surface.include(curve),
                                            Curve::ParameterCurve(_) => false,
                                            Curve::IntersectionCurve(_) => false,
                                        }
                                    }
                                    Curve::ParameterCurve(_) => unreachable!(),
                                }
                            }
                            Curve::ParameterCurve(_) => false,
                            Curve::IntersectionCurve(_) => false,
                        }
                    }
                    Curve::ParameterCurve(_) => false,
                },
            }
        }
    }
}

impl IncludeCurve<Curve> for Plane {
    fn include(&self, curve: &Curve) -> bool {
        curve.lift_up().control_points().iter().all(|v| {
            let p = v.to_point();
            self.search_parameter(p, None, 1).is_some()
        })
    }
}

impl ToSameGeometry<Surface> for Plane {
    fn to_same_geometry(&self) -> Surface { (*self).into() }
}

// Lossless extraction of a plane, used by planar-only algorithms such as
// `monstertruck_solid::shell` and `monstertruck_solid::draft`.
impl TryFrom<Surface> for Plane {
    type Error = ();
    fn try_from(surface: Surface) -> std::result::Result<Self, ()> {
        match surface {
            Surface::Plane(plane) => Ok(plane),
            _ => Err(()),
        }
    }
}

impl ToSameGeometry<Surface> for RevolutionSurface<Curve> {
    fn to_same_geometry(&self) -> Surface {
        Surface::RevolutionSurface(Processor::new(self.clone()))
    }
}

impl SearchNearestParameter<SurfaceParameter> for Surface {
    type Point = Point3;
    fn search_nearest_parameter<H: Into<SearchParameterHint2D>>(
        &self,
        point: Point3,
        hint: H,
        trials: usize,
    ) -> Option<(f64, f64)> {
        match self {
            Surface::Plane(plane) => plane.search_nearest_parameter(point, hint, trials),
            Surface::BsplineSurface(bspsurface) => {
                bspsurface.search_nearest_parameter(point, hint, trials)
            }
            Surface::NurbsSurface(surface) => surface.search_nearest_parameter(point, hint, trials),
            Surface::TsplineSurface(surface) => {
                surface.search_nearest_parameter(point, hint, trials)
            }
            Surface::RevolutionSurface(rotted) => {
                let hint = match hint.into() {
                    SearchParameterHint2D::Parameter(hint0, hint1) => (hint0, hint1),
                    SearchParameterHint2D::Range(x, y) => {
                        algo::surface::presearch(rotted, point, (x, y), 100)
                    }
                    SearchParameterHint2D::None => {
                        algo::surface::presearch(rotted, point, rotted.range_tuple(), 100)
                    }
                };
                algo::surface::search_nearest_parameter(rotted, point, hint, trials).or_else(|| {
                    let candidate = rotted.evaluate(hint.0, hint.1);
                    if candidate.near(&point) {
                        Some(hint)
                    } else {
                        None
                    }
                })
            }
        }
    }
}

impl TryIntoBsplineSurface for Surface {
    fn try_into_bspline_surface(&self) -> Option<BsplineSurface<Point3>> {
        match self {
            Surface::Plane(p) => p.try_into_bspline_surface(),
            Surface::BsplineSurface(b) => b.try_into_bspline_surface(),
            Surface::NurbsSurface(n) => n.try_into_bspline_surface(),
            Surface::RevolutionSurface(_) => None,
            Surface::TsplineSurface(_) => None,
        }
    }
}

impl TryIntoHomogeneousBsplineSurface for Surface {
    fn try_into_homogeneous_bspline_surface(&self) -> Option<BsplineSurface<Vector4>> {
        match self {
            Surface::Plane(p) => p.try_into_homogeneous_bspline_surface(),
            Surface::BsplineSurface(b) => b.try_into_homogeneous_bspline_surface(),
            Surface::NurbsSurface(n) => n.try_into_homogeneous_bspline_surface(),
            Surface::RevolutionSurface(r) => r.try_into_homogeneous_bspline_surface(),
            Surface::TsplineSurface(_) => None,
        }
    }
}

impl SupportsExactPatchDomains for Surface {
    fn supports_exact_patch_domains(&self) -> bool {
        matches!(self, Surface::BsplineSurface(_) | Surface::NurbsSurface(_))
    }
}

impl ToSameGeometry<Surface> for HomotopySurface<Curve, Curve> {
    fn to_same_geometry(&self) -> Surface {
        let curve0 = self.first_curve().clone().lift_up();
        let curve1 = self.second_curve().clone().lift_up();
        NurbsSurface::new(BsplineSurface::homotopy(curve0, curve1)).into()
    }
}

impl ToSameGeometry<Surface> for ExtrusionSurface<Curve, Vector3> {
    fn to_same_geometry(&self) -> Surface {
        let (curve0, vector) = (self.entity_curve(), self.extruding_vector());
        let trsl = Matrix4::from_translation(vector);
        let curve1 = self.entity_curve().transformed(trsl);
        match (curve0, curve1) {
            (Curve::Line(line), Curve::Line(_)) => {
                Plane::new(line.0, line.1, line.0 + vector).into()
            }
            (Curve::BsplineCurve(curve0), Curve::BsplineCurve(curve1)) => {
                BsplineSurface::homotopy(curve0.clone(), curve1.clone()).into()
            }
            (Curve::NurbsCurve(curve0), Curve::NurbsCurve(curve1)) => {
                NurbsSurface::new(BsplineSurface::homotopy(
                    curve0.non_rationalized().clone(),
                    curve1.non_rationalized().clone(),
                ))
                .into()
            }
            (Curve::IntersectionCurve(_), Curve::IntersectionCurve(_)) => unimplemented!(),
            _ => unreachable!(),
        }
    }
}

// ---------------------------------------------------------------------------
// Deterministic content hashing for modeling enums.
// ---------------------------------------------------------------------------

impl monstertruck_core::DeterministicContentHash for Conic2D {
    fn content_hash<H: std::hash::Hasher>(&self, state: &mut H) {
        match self {
            Self::Ellipse(v) => {
                state.write_u8(0);
                v.content_hash(state);
            }
            Self::Hyperbola(v) => {
                state.write_u8(1);
                v.content_hash(state);
            }
            Self::Parabola(v) => {
                state.write_u8(2);
                v.content_hash(state);
            }
        }
    }
}

impl monstertruck_core::DeterministicContentHash for Curve2D {
    fn content_hash<H: std::hash::Hasher>(&self, state: &mut H) {
        match self {
            Self::Line(v) => {
                state.write_u8(0);
                v.content_hash(state);
            }
            Self::Polyline(v) => {
                state.write_u8(1);
                v.content_hash(state);
            }
            Self::Conic(v) => {
                state.write_u8(2);
                v.content_hash(state);
            }
            Self::BsplineCurve(v) => {
                state.write_u8(3);
                v.content_hash(state);
            }
            Self::NurbsCurve(v) => {
                state.write_u8(4);
                v.content_hash(state);
            }
        }
    }
}

impl monstertruck_core::DeterministicContentHash for Curve {
    fn content_hash<H: std::hash::Hasher>(&self, state: &mut H) {
        match self {
            Self::Line(v) => {
                state.write_u8(0);
                v.content_hash(state);
            }
            Self::BsplineCurve(v) => {
                state.write_u8(1);
                v.content_hash(state);
            }
            Self::NurbsCurve(v) => {
                state.write_u8(2);
                v.content_hash(state);
            }
            Self::ParameterCurve(v) => {
                state.write_u8(3);
                v.content_hash(state);
            }
            Self::IntersectionCurve(v) => {
                state.write_u8(4);
                v.content_hash(state);
            }
        }
    }
}

impl monstertruck_core::DeterministicContentHash for Surface {
    fn content_hash<H: std::hash::Hasher>(&self, state: &mut H) {
        match self {
            Self::Plane(v) => {
                state.write_u8(0);
                v.content_hash(state);
            }
            Self::BsplineSurface(v) => {
                state.write_u8(1);
                v.content_hash(state);
            }
            Self::NurbsSurface(v) => {
                state.write_u8(2);
                v.content_hash(state);
            }
            Self::RevolutionSurface(v) => {
                state.write_u8(3);
                v.content_hash(state);
            }
            Self::TsplineSurface(v) => {
                state.write_u8(4);
                v.content_hash(state);
            }
        }
    }
}
