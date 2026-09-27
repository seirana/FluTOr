function coupling = computeCouplings(model, options)
%COMPUTECOUPLINGS Compute QFCA/F2C2 consensus full-coupling evidence.
%
%   OPTIONS:
%     tolerance     QFCA numerical tolerance (default 1e-10)
%     f2c2Solver    solver name passed to F2C2 (default "glpk")
%     f2c2Matrix    optional precomputed F2C2 full-coupling matrix
%
%   A precomputed matrix is useful for reproducible environments where the
%   external F2C2 implementation is unavailable in CI.

    if nargin < 2 || isempty(options)
        options = struct();
    end
    if ~isfield(options, "tolerance")
        options.tolerance = 1e-10;
    end
    if ~isfield(options, "f2c2Solver")
        options.f2c2Solver = "glpk";
    end
    if ~isfield(options, "f2c2Matrix")
        options.f2c2Matrix = [];
    end

    model = flutor.validateModel(model);
    [qfcaMatrix, qfcaDiagnostics] = ...
        flutor.computeQfcaFullCoupling(model, options.tolerance);

    if isempty(options.f2c2Matrix)
        [f2c2Matrix, f2c2Diagnostics] = ...
            flutor.computeF2C2FullCoupling(model, options.f2c2Solver);
    else
        f2c2Matrix = logical(options.f2c2Matrix);
        f2c2Diagnostics = struct( ...
            "solver", "precomputed", ...
            "nFullCoupledPairs", nnz(triu(f2c2Matrix, 1)));
    end

    coupling = flutor.combineFullCouplingEvidence( ...
        qfcaMatrix, ...
        f2c2Matrix);
    coupling.diagnostics = struct( ...
        "qfca", qfcaDiagnostics, ...
        "f2c2", f2c2Diagnostics);
end
