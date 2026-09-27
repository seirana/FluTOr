# Modernization notes

The publication-era implementation is preserved under `legacy/`. The maintained package is under `src/+flutor/`.

## Main engineering changes

The modernization replaces script-level state with explicit function APIs.

Historical pattern:

```text
workspace variables
  -> load model
  -> mutate bounds
  -> linprog / COBRA calls
  -> coupling kernels
  -> MILP loop
  -> save files during enumeration
```

Maintained pattern:

```text
validated model
  -> explicit bound-change table
  -> checked preprocessing
  -> explicit coupling evidence
  -> checked MILP
  -> structured result
  -> separate persistence
```

## Concrete defects addressed

### Undeclared workspace variables

Several historical experiment scripts use variables such as `adr`, `adrs`, or `f_name` without defining them in the script itself.

The maintained API passes paths and run names as function arguments.

### Solver failure treated as zero flux

Historical FBA/FVA code assigned zero when `linprog` did not return an optimal status.

The maintained implementation raises `flutor:OptimizationFailed` instead. Infeasibility or solver failure is not equivalent to a biologically blocked reaction.

### Reduced reaction-count mutation in QFCA

The historical modified QFCA decremented its reduced reaction-count variable as soon as a metabolite row contained two nonzeros, before confirming that the reactions were actually merged.

The maintained full-coupling kernel updates dimensions only after a successful merge.

### External F2C2 dependency

The historical code invokes `F2C2` but the implementation is not stored in this repository.

The maintained adapter treats F2C2 as an explicit external dependency and validates its output dimensions.

### Output during enumeration

The historical trade-off seeker saves the MAT result repeatedly from inside the optimization loop.

The maintained enumerator returns structured results; `flutor.writeResults` performs persistence after computation.

## Scientific boundary

The original source remains available for publication-era reproduction. If maintained and historical results differ, compare identical:

- model snapshots;
- environmental constraints;
- biomass target;
- solver versions;
- F2C2 implementation;
- numerical tolerances.

Do not assume a difference is caused by code quality changes before controlling these factors.
