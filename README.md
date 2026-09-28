# FluTOr — Relative Flux Trade-Offs in Metabolic Networks

FluTOr is a MATLAB implementation of a constraint-based method for identifying relative flux trade-offs with respect to an optimized metabolic task such as growth.

The method is associated with:

> Hashemi, S., Razaghi-Moghadam, Z., Laitinen, R. A. E., & Nikoloski, Z. (2022). **Relative flux trade-offs and optimization of metabolic network functionalities.** *Computational and Structural Biotechnology Journal*, 20, 3963–3971. https://doi.org/10.1016/j.csbj.2022.07.038

## Maintained versus historical code

The repository now separates the original publication-era scripts from the maintained implementation:

```text
src/+flutor/       maintained package
legacy/            original scripts and functions
```

New analyses should use the maintained package.

## What FluTOr studies

FluTOr seeks relations in which the flux of a target reaction, typically biomass, is represented by a weighted combination of other reaction fluxes under a constrained metabolic network.

The maintained workflow is:

```text
validated metabolic model
      |
      v
explicit reaction-bound changes
      |
      v
checked FBA/FVA preprocessing
      |
      v
irreversible representation
      |
      v
QFCA full-coupling evidence
      +
F2C2 full-coupling evidence
      |
      v
consensus coupling groups
      |
      v
relative-trade-off MILP
      |
      v
structured outputs + provenance
```

## Why the modernization was necessary

The original implementation contains several characteristics of exploratory research code:

- top-level scripts depended on undeclared workspace variables such as output paths;
- solver configuration and model manipulation were mixed with analysis logic;
- FBA/FVA code converted failed optimizations into zero flux bounds;
- numerical equality relied on decimal rounding;
- F2C2 was called as an undeclared external dependency;
- the modified QFCA kernel contained mutable loop state that was difficult to audit;
- one reaction-count decrement occurred before confirming that a merge happened;
- output files were written from inside the enumeration loop;
- there were no automated tests or CI.

The maintained package addresses those engineering problems while preserving the publication-era source under `legacy/`.

## Maintained package

```text
src/+flutor/
├── validateModel.m
├── applyReactionBounds.m
├── computeFluxRanges.m
├── computeQfcaFullCoupling.m
├── computeF2C2FullCoupling.m
├── combineFullCouplingEvidence.m
├── computeCouplings.m
├── expandCoupledAlternatives.m
├── enumerateRelativeTradeoffs.m
├── runAnalysis.m
└── writeResults.m
```

The code follows a function-oriented design with explicit inputs/outputs, namespaced errors, numerical tolerances, checked solver status, and separation of computation from persistence.

## Dependencies

Full analysis requires:

- MATLAB;
- Optimization Toolbox (`linprog`, `intlinprog`);
- COBRA Toolbox (`readCbModel`, `convertToIrreversible`, model-removal helpers);
- the F2C2 implementation used by the original workflow;
- an LP solver supported by that F2C2 installation.

The repository does not vendor F2C2. `flutor.computeF2C2FullCoupling` detects its absence and fails with an explicit dependency error.

For fully reproducible runs, record the exact F2C2 source/version and solver version.

## Setup

```matlab
addpath("src")
```

Initialize COBRA Toolbox and F2C2 separately.

## Model validation

`flutor.validateModel` requires:

```text
S
rxns
mets
lb
ub
c
```

It validates dimensions, finite numerical values, and reaction bounds. Stable reaction/metabolite indices are created if missing.

## Replacing hard-coded condition scripts

Instead of embedding model-specific changes in a script, create a table:

```matlab
changes = table( ...
    ["EX_glc__D_e"; "BIOMASS_Ec_iJO1366_WT_53p95M"], ...
    [-10; 0.8], ...
    [1000; 1000], ...
    'VariableNames', ...
    {'ReactionID', 'LowerBound', 'UpperBound'});

model = flutor.applyReactionBounds(model, changes);
```

This makes the biological condition auditable and reusable.

## FBA/FVA preprocessing

`flutor.computeFluxRanges`:

- maximizes the target/biomass reaction when a biomass fraction is requested;
- applies optional additional linear inequalities;
- computes per-reaction minima and maxima with `linprog`;
- checks every optimization exit flag;
- does **not** silently replace solver failure with zero;
- converts the model to an irreversible representation through COBRA Toolbox;
- removes blocked reactions;
- records preprocessing diagnostics.

Example:

```matlab
options = struct( ...
    "biomassFraction", 0.90, ...
    "tolerance", 1e-5);

[processedModel, diagnostics] = ...
    flutor.computeFluxRanges( ...
        model, ...
        "BIOMASS_Ec_iJO1366_WT_53p95M", ...
        options);
```

## Coupling analysis

FluTOr historically compared full-coupling calls from modified QFCA and F2C2.

The maintained code makes that consensus explicit:

```matlab
coupling = flutor.computeCouplings(processedModel);
```

The output contains:

- QFCA full-coupling matrix;
- F2C2 full-coupling matrix;
- intersection/consensus matrix;
- deterministic coupling groups;
- diagnostics.

A precomputed F2C2 matrix can be supplied through `couplingOptions.f2c2Matrix` for controlled or offline environments.

## Relative trade-off enumeration

`flutor.enumerateRelativeTradeoffs` constructs the MILP using sparse matrices and:

- validates the target reaction;
- excludes the biomass full-coupling group from candidate variable reactions;
- prevents simultaneous selection of fully-coupled alternatives;
- checks `intlinprog` termination;
- distinguishes infeasibility from solver failure;
- expands coupling-equivalent supports deterministically;
- removes duplicate supports;
- returns a table instead of saving during every iteration.

## End-to-end API

```matlab
addpath("src")

model = readCbModel("models/Ecoli_iJO1366.mat");

options = struct();
options.fluxRangeOptions = struct( ...
    "biomassFraction", 0.90);

output = flutor.runAnalysis( ...
    model, ...
    "BIOMASS_Ec_iJO1366_WT_53p95M", ...
    options);
```

## File-driven runner

```matlab
addpath("scripts")

output = runFluTOr( ...
    "models/Ecoli_iJO1366.mat", ...
    "BIOMASS_Ec_iJO1366_WT_53p95M", ...
    "artifacts/ecoli");
```

An optional CSV of reaction-bound changes can be provided as the fourth argument.

The runner records SHA-256 checksums for the model and optional bounds file.

## Tests

```matlab
addpath("src")
results = runtests("tests");
assertSuccess(results)
```

Solver-independent tests cover:

- model validation;
- reaction-ID bound changes;
- QFCA/F2C2 consensus behavior;
- deterministic coupling-alternative expansion;
- a synthetic QFCA fully-coupled pair.

CI does not require COBRA or F2C2 because external scientific dependencies should not make core validation logic untestable.

## Repository layout

```text
.
├── src/+flutor/        maintained package
├── tests/              unit tests
├── scripts/            reproducible runner
├── legacy/             publication-era implementation
├── models/             model snapshots and condition tables
├── Results/            historical publication output
├── README.md
├── MIGRATION.md
├── REPRODUCIBILITY.md
└── CITATION.cff
```

## Scientific limitations

Results are conditional on the metabolic reconstruction, constraints, target reaction, coupling method, numerical tolerances, solver behavior, and environmental condition.

Computationally identified relative trade-offs and overexpression candidates are hypotheses derived from the model. They are not experimental validation.

## License

The repository already contains an MIT license. See [LICENSE](LICENSE).
