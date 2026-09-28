# Reproducibility

A FluTOr experiment should preserve:

- repository Git commit;
- MATLAB release;
- Optimization Toolbox version through the MATLAB release;
- COBRA Toolbox version/commit;
- F2C2 source/version;
- LP solver and version;
- model source/version/checksum;
- reaction-bound condition file;
- biomass/target reaction;
- biomass fraction;
- additional inequality constraints;
- numerical tolerances;
- MILP coefficient/Big-M limits;
- MILP time limit;
- generated trade-off table and analysis MAT file.

## External F2C2

F2C2 is not vendored by this repository.

For publication-scale work, archive either:

1. the exact F2C2 implementation used, subject to its license; or
2. the computed F2C2 full-coupling matrix with provenance and checksum.

The maintained `flutor.computeCouplings` accepts a precomputed F2C2 matrix.

## Model-condition provenance

Avoid relying on row numbers from spreadsheets unless the exact model snapshot is also frozen.

Prefer reaction IDs in condition tables whenever possible.

## Solver failures

A solver failure must remain visible in the experiment record. The maintained code does not reinterpret a non-optimal solver status as zero flux.

## Historical comparison

The original scripts are under `legacy/`. Use them only with their original dependency stack when reproducing publication-era output.
