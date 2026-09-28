# Contributing

## Maintained code

New development belongs under:

```text
src/+flutor/
```

Publication-era scripts under `legacy/` are retained for provenance and should not be rewritten except for archival/documentation maintenance.

## Development workflow

1. create a branch;
2. keep scientific and engineering changes reviewable;
3. add or update MATLAB unit tests;
4. run the test suite;
5. update README/MIGRATION/REPRODUCIBILITY documentation when behavior changes;
6. open a pull request.

## Tests

```matlab
addpath("src")
results = runtests("tests");
assertSuccess(results)
```

Prefer small synthetic models for unit tests so CI does not depend on COBRA Toolbox, F2C2, or a configured external solver unless the behavior under test genuinely requires them.

## Scientific changes

Any change that may alter FluTOr results should document:

- previous behavior;
- new behavior;
- reason for the change;
- expected effect on publication-era results;
- model/solver assumptions;
- validation performed.

## Style

Prefer:

- package-qualified public functions;
- descriptive names;
- explicit options structs;
- validated model dimensions and bounds;
- numerical tolerances instead of equality-by-rounding;
- sparse matrices for optimization formulations;
- checked solver exit flags;
- deterministic ordering;
- structured return values;
- file I/O only at orchestration boundaries.
