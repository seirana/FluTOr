# Legacy publication implementation

This directory contains the original FluTOr MATLAB scripts and functions.

They are preserved for provenance and publication-era comparison, not as the maintained API.

Known engineering limitations include:

- undeclared workspace variables in experiment scripts;
- model- and path-specific assumptions;
- direct mutation of global solver/path state;
- solver failures converted to zero flux bounds;
- external F2C2 dependency not vendored in the repository;
- output persistence inside optimization loops;
- limited automated validation and no tests.

Use `src/+flutor/` for new work.

See [../MIGRATION.md](../MIGRATION.md).
