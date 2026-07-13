# Design note: control references in algebraic expressions

Status: design, not yet implemented. Records the fix agreed 2026-07-13.

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

### Secondary: the final-time control

A control decision exists for elements `0..N-1`; there is none at the final
node `fe[N]` (nothing evolves after the final time). So an **algebraic**
reference to the final-time control should error rather than be silently
rewritten. For piecewise-constant that variable does not survive
`parameterize` anyway; the value of the fix is turning a masked mistake
into a clear message. This is a safety net around the correctness fix, not
the fix itself, and it does not touch the differential equations (which
legitimately reference the final node — e.g. piecewise-linear's last-
segment interpolation).

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
