# Design note: control references in algebraic expressions

Status: implemented and tested 2026-07-13 (the substitution split in
`parameterize.py`; `tests/test_algebraic_refs.py`). This note records the
design and the profile x scheme coverage analysis behind it.

## Summary

The bug is a correctness bug, and it is confined to **piecewise-constant**
controls in **algebraic** expressions (objectives, tracking constraints,
terminal costs, path constraints). `cvp.parameterize` rewrites control
references the same way in every expression. Inside the differential
equations that is correct. Inside an algebraic expression it is wrong: a
reference to element k's control, `v[fe[k]]`, is rewritten to element k-1's
control, `v[fe[k-1]]`. The objective then penalizes the wrong element (and
double-weights element 0). It is a silent wrong answer, not an error.

The fix: rewrite control references **only inside the differential
equations**. In algebraic expressions, a control reference is the literal
element control, `v[fe[k]] = element k's control`, with no shift. The
dynamics and pyomo.dae are untouched.

Erroring on a reference to the (nonexistent) final-time control `v[fe[N]]`
falls out of the same change and is a useful safety net, but it is
secondary. The point is that the controls the objective references are the
ones the user meant.

## What is already correct and does not change

- **The ODEs (all profiles).** Verified: the differential equation at a
  finite-element boundary `fe[k]` uses `v[fe[k-1]]` (element k-1's control,
  because `fe[k]` is the right-hand collocation node of element k-1);
  interior collocation points of element k use `v[fe[k]]`. The equation at
  the final node `fe[N]` uses `v[fe[N-1]]`. All correct. No change.
- **piecewise_linear (all schemes).** Its algebraic references do not shift
  (a reference at `fe[k]` maps to its own knot `v[fe[k]]`), and its
  dynamics correctly use `v[fe[N]]` as the interpolation endpoint the last
  segment ramps toward (`v(t) = (1-w)·v[fe[N-1]] + w·v[fe[N]]` on the last
  element, verified). `v[fe[N]]` exists and is used by the dynamics, and is
  not a control the user would reference. Leave it alone.
- **reduced_collocation.** Structurally different (free values at
  collocation points, not element boundaries; interior points Lagrange-
  interpolated). Out of scope for this fix.

## The bug (piecewise_constant, algebraic expressions)

For piecewise-constant control there are N controls, one per element, and
after `parameterize` the surviving variable is indexed at the element
starts `fe[0..N-1]`. `v[fe[k]]` is element k's control.

The transform builds one substitution map and applies it to every
constraint and every objective (`parameterize.py`, the loops over
`Constraint` and `Objective` that call `replace_expressions(..., submap)`).
That map sends a control reference at time t to the control the ODE at t
uses. At a boundary `fe[k]` that is element k-1's control, `v[fe[k-1]]`.

- In a differential equation that is exactly right: the equation at `fe[k]`
  is element k-1's right collocation node.
- In an algebraic expression it is wrong. A user writing `v[fe[k]]` in an
  objective means element k's control. The transform rewrites it to element
  k-1's. Seen in the Hicks example (`examples/hicks_cvp.ipynb`): a tracking
  objective `sum((v[i] - uss)**2 for i in m.i)` becomes
  `v[fe[0]], v[fe[0]], v[fe[1]], ..., v[fe[N-1]]` — element 0 penalized
  twice, every other element shifted down by one. Wrong objective, no error.

The reference convention is the AMPL Hicks model
(`~/Dropbox/CMU Research/Hicks/hicks.mod`): controls are element-indexed
(`u1{i in fe}`, `fe = 0..N-1`), tracking references `u1[i]` per element,
and no node is shared, so the shift cannot arise there.

## The fix

Split the substitution by expression kind, distinguished by whether the
expression references a `DerivativeVar`:

- **Differential equations** (reference a `DerivativeVar`): unchanged. Keep
  the current map. Correct for every profile and scheme.
- **Algebraic expressions** (no `DerivativeVar`): a control reference maps
  to the literal element control. For piecewise-constant, `v[fe[k]]` stays
  `v[fe[k]]` (element k's control), not `v[fe[k-1]]`. For piecewise-linear
  and reduced_collocation the algebraic map is unchanged (no shift exists
  to correct).

After the fix a tracking objective summed over the element starts
`fe[0..N-1]` penalizes each element's control once, at its own index,
matching the AMPL `sum{i in fe} ... u1[i]`.

**Detection gotcha.** "References a `DerivativeVar`" must be tested by
walking `identify_variables(expr)` and checking membership in the set of
DerivativeVar data ids (or `isinstance(v.parent_component(), DerivativeVar)`).
`identify_components(expr, [DerivativeVar])` returns nothing, a false
negative, so it cannot be used. After `dae.collocation` or
`dae.finite_difference` the user's ODE constraint still carries its
`DerivativeVar` (e.g. `zdot[t] == u[t] - z[t]`), and pyomo.dae's own
`disc_eq` rows carry it too; both are differential and both take the
differential map (the `disc_eq` reference no controls, so it is a no-op
there).

### The final-time control: an error, not a mapping

For piecewise-constant there is no control at the final node `fe[N]`
(nothing evolves past the final time), and the variable does not survive
`parameterize`. An algebraic reference to it cannot map to a real element
control. The only non-erroring option would be to map it to the last
element `v[fe[N-1]]`, but that silently double-counts the last element's
control, which is exactly the silent wrong answer this fix removes. So the
final node must error, with a message naming the control and pointing at
the elements. This is part of correctness, not a separate safety feature.
It does not touch the differential equations, which legitimately reference
the final node (piecewise-linear's last-segment interpolation, and every
profile's final collocation ODE).

## The real models (checked)

Both examples route the objective through a scalar `Var` `m.track` defined
by an algebraic `Constraint` `tracking_def` (`m.track == sum(...)`), the
AMPL `tracking`/`trackingdef` pattern, not a named `Expression`. So:

- The controls live in `tracking_def`, a `DerivativeVar`-free constraint,
  correctly classified algebraic. The fix reaches them.
- **quad-tank** tracks states only, no controls, so it is unaffected.
- **hicks** tracks states and controls and sums over all of `m.i`, the
  boundaries `{0..N}`, so `tracking_def` references `v1[fe[N]]`, `v2[fe[N]]`,
  the final-node controls. Under this fix that constraint will error until
  the tracking stops penalizing the terminal control. The AMPL-consistent
  form sums the control terms over the elements `0..N-1` and handles the
  terminal state separately (AMPL's `termcost`, which references no
  controls). hicks_cvp.ipynb is untracked WIP and is not edited here; the
  change is on the user's side.

No control-referencing named `Expression`s appear in the real models. The
implementation still classifies named `Expression`s by their own
`DerivativeVar` content; a control-referencing `Expression` shared between a
differential and an algebraic root would resolve by that content and is a
documented limitation, not a case the real models hit.

## Scope and non-goals

- No change to pyomo.dae, the discretization, or the differential
  equations. The dynamics are already correct.
- No change to the control's declaration surface (`declare_profile`) or the
  global-time indexing. Controls stay `Var(ContinuousSet)`; this is only
  about how references are rewritten.
- piecewise_linear and reduced_collocation algebraic references are
  unchanged.

## Tests to add

- piecewise-constant, objective referencing controls at the element starts:
  after `parameterize`, each element control appears once at its own index
  (no shift, no double-count). This is the correctness test.
- Regression: the differential path still rewrites the boundary node to
  `v[fe[k-1]]` and the final node to `v[fe[N-1]]` (differential handling
  unchanged).
- piecewise-linear objective referencing controls at all boundaries is
  unchanged by the fix (no shift to correct), and its dynamics still use
  `v[fe[N]]`.
- Secondary: a piecewise-constant algebraic reference to the final-time
  control raises a clear error.
