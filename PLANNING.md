# Kernel-side planning (fork workspace docs — never merged upstream)

Companion to CONVENTIONS.md. Downstream state, runtime plans, and the status audit live
in the openshape repo (`notes/` — see its `notes/handoff/`). This file holds the
kernel-side project briefs. Existing artifacts are referenced, not duplicated: branch
inventory in CONVENTIONS.md; per-op capability matrix in openshape
`notes/compatibility.md`; delivered PR/issue bodies in the openshape conversation log
and the branches themselves.

## Brief 1 — Tangent/coincident-face booleans (real solution)

Status quo: `tr/fix/coincident-plane-booleans` fixes the two sharpest edges (panic →
`None`; coplanar seg/triangle NaN) and should be upstreamed FIRST as-is — it lands
independently of any design work. The remaining gap (touching-face union still `None`)
is a design problem, in truck-shapeops' terms: "faces tangent to each other are not yet
supported."

Design direction (from the diagnosis on that branch): the failure chain is
(1) coincident-plane face pairs produce no transversal intersection curve, so
`divide_face` never partitions the shared region; (2) spurious near-zero loops from
perpendicular-edge intersections get spliced into the shared face's loop list;
(3) `PolylineCurve::include`'s ray-cast returns false for on-vertex points.
A real fix needs, in order: (a) detect coincident/near-coincident (within tolerance)
face pairs up front; (b) for each pair, run 2D polygon-overlay classification (clip the
two faces' loop sets against each other in the shared plane — union/intersection/
difference regions labeled per boolean op); (c) replace the pair with the classified
regions before the normal transversal pipeline runs; (d) make containment tests
boundary-robust (on-edge/on-vertex → deterministic tie-break, not `unwrap_or(false)`).
The 2D overlay is standard computational geometry (e.g. Greiner–Hormann or
Martinez–Rueda) and could vendor a small robust implementation or be written against
truck's own 2D polyline types. Size: L. Do it against BOTH upstreams? Write for truck
(`tr/feat/coplanar-boolean-classification`) since truck-shapeops is the reference
implementation monstertruck reverted TO; port to monstertruck after review feedback.

## Brief 2 — Curved-face shell/thicken/draft

Investigated on `mt/feat/planar-shared-and-curved` (shipped the planar.rs refactor
only; see that branch's report in the openshape log). Precise blocker recorded there:
`ThickenableSurface/Curve` generics hardcode lossless `Plane`/`Line` round-trips, and
plane∩cylinder corners are quadratic multi-root problems. The real project:
1. Introduce solid-local sum types `ShellSurface {Plane, Cylinder(...)}` /
   `ShellCurve {Line, Arc(...)}` (cylinders = `RevolutionSurface<Line>`; arcs =
   `Processor<TrimmedCurve<UnitCircle>>` — both already in monstertruck-geometry).
2. Re-generic the engine over those sum types with total `From` back-conversions.
3. Corner solver: extend `planar.rs` with mixed constraint sets — plane/plane/cylinder
   and plane/cylinder/cylinder corner points (quadratic; disambiguate roots by
   proximity to the original vertex), cylinder-parallel edge offsets, arc offsets.
4. Validation: cylinder hollow (analytic pi·(r²h − (r−t)²(h−t')) cases), rounded-corner
   box, half-pipe; reject non-concentric cylinder-cylinder corners initially.
Size: L (shell/thicken) + M (extending draft similarly). Base on `patches/monstertruck`
until shell+draft land upstream.

## Brief 3 — Sketching (2D) in the kernel

Facts: upstream truck has `truck-drafting` (2D drafting crate, growing: line/arc
connectors). monstertruck deferred it deliberately — TRUCK-PARITY.md: "not
public-API quality yet (panicking wrappers, ambiguous arc-constraint names,
scale-dependent arc-length integration, broad prelude re-exports)"; its cleanup
checklist is in truck-sync.md. "No sketching in the kernel" does NOT mean no 2D basis
for solids: profiles are built directly as wires (vertices/lines/arcs/Béziers via
builder::*), which is exactly what openshape's runtime does (polygon profiles today).
What a sketch layer adds: arcs/circles/splines in profiles with region extraction
(closed-loop detection, face assembly from mixed wire sets) and, separately,
CONSTRAINT SOLVING — which is NOT a kernel concern in Onshape either (see openshape
notes/handoff for the FeatureScript-side answer). Plan: (a) port/clean truck-drafting
into monstertruck-sketch per the existing checklist (M, mostly mechanical + API
hygiene); (b) region extraction from wire soup (M); constraint solving stays out of
the kernel.

## Brief 4 — other kernel gaps

See the audit table in openshape `notes/handoff/status-audit.md` (patterns/mirror =
transform ops that mostly exist; helix; surface ops; NURBS-surface creation from
FeatureScript-shaped inputs). Triage there.
