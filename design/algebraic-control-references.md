# Design note: control references in algebraic expressions

Status: design, not yet implemented. Records the fix agreed 2026-07-13.

## Summary

`cvp.parameterize` rewrites every control reference in the model uniformly.
That is correct inside the differential (ODE/collocation) equations and
wrong inside algebraic expressions (objectives, tracking constraints,
terminal costs, path constraints). The fix is to split the substitution
by expression kind: keep the current rewrite for differential equations,
and in algebraic expressions treat a control reference as the literal
control of the element it names, erroring when it names the final node.

pyomo.dae and the discretization are untouched. The dynamics are already
correct.

## Background: how the control is indexed

A control is a `Var` over the discretized `ContinuousSet`, so it exists at
every time node: the finite-element boundaries `fe[0..N]` and the interior
collocation points. Piecewise-constant control has N free values, one per
element. After `parameterize` the surviving control variable is indexed at
the element starts `fe[0..N-1]`:

- `v[fe[k]]` is element k's control (the value held across element k).
- `v[fe[N]]` does not exist. There are N elements, hence N controls.

The control that the ODE at a given time uses:

- At an interior collocation point of element k: element k's control,
  `v[fe[k]]`.
- At a finite-element boundary `fe[k]`: element k-1's control, `v[fe[k-1]]`,
  because `fe[k]` is the right-hand collocation node of element k-1 (Radau
  places a node at each element's right endpoint). The ODE at the final
  node `fe[N]` therefore uses `v[fe[N-1]]`.

So `v[fe[N]]` is never a legitimate quantity: the differential equation at
the final node references `v[fe[N-1]]`, and there is no element N.

## The problem

`parameterize` replaces control references the same way in every
expression (see `parameterize.py`, the loops over `Constraint` and
`Objective` that call `replace_expressions(..., submap)`). The `submap`
sends each control reference to the control the ODE at that time uses. That
is exactly right for the differential equations.

It is wrong for algebraic expressions, because a user writing a control in
an objective means "element k's control", not "the control the ODE at
`fe[k]` uses". Two symptoms, both seen in the Hicks example
(`examples/hicks_cvp.ipynb`):

1. A control-tracking objective written as
   `sum((v[i] - uss)**2 for i in m.i)` references controls at the
   finite-element boundaries. Each `v[fe[k]]` is rewritten to `v[fe[k-1]]`
   (the ODE-at-`fe[k]` value), so the objective penalizes the wrong element
   and double-counts `v[fe[0]]`. The user meant element k's control.
2. The objective also references `v[fe[N]]`. Because the reference was
   written before `parameterize` (when the control still existed at every
   node), it is silently rewritten to `v[fe[N-1]]` instead of erroring. The
   mistake — referencing a control that does not exist — is masked.

The reference model is `~/Dropbox/CMU Research/Hicks/hicks.mod` (AMPL),
where controls are indexed by element (`u1{i in fe}`, `fe = 0..N-1`):
tracking references `u1[i]` per element, `u1[N]` cannot be named, and the
dynamics reference `u2[i]` by element. No node is shared, so none of this
arises.

## The fix

Split the substitution pass by expression kind, distinguished by whether
the expression references a `DerivativeVar`:

- **Differential equations** (reference a `DerivativeVar`): unchanged.
  Keep the current rewrite. The equation at `fe[N]` continues to come out
  as `v[fe[N-1]]`, which is correct.
- **Algebraic expressions** (no `DerivativeVar`: objectives, tracking and
  terminal constraints, path constraints): a control reference maps to the
  literal control of the element it names — `v[fe[k]]` stays `v[fe[k]]`
  (element k's control), not `v[fe[k-1]]`. A reference at the final node
  `fe[N]` raises a clear error naming the final-time control, because no
  element owns it.

After the fix, a control-tracking objective summed over the element starts
`fe[0..N-1]` penalizes each element's control once, matching the AMPL
`sum{i in fe} ... u1[i]`, and any stray reference to the final-time control
fails loudly instead of being rewritten.

## Scope and non-goals

- No change to pyomo.dae, the discretization, or the differential
  equations. The dynamics are already correct.
- No change to the control's declaration surface (`declare_profile`) or to
  the global-time indexing. Controls remain `Var(ContinuousSet)`; this is
  purely about how references are rewritten.
- Interior-collocation-point control references inside algebraic
  expressions (unusual) map to the containing element's control; only the
  final node is an error.

## Tests to add

- A model whose objective references controls at the element starts:
  after `parameterize`, each element control appears once, at its own
  index (no shift, no double-count).
- A model whose objective references the final-time control `v[fe[N]]`:
  `parameterize` raises, naming the final-time control.
- A differential equation referencing `v[fe[N]]` (the standard ODE written
  over the whole set) still parameterizes cleanly to `v[fe[N-1]]`
  (regression: the differential path is unchanged).
