function [fullCoupled, diagnostics] = computeF2C2FullCoupling(model, solverName)
%COMPUTEF2C2FULLCOUPLING Run an installed F2C2 implementation safely.
%
%   F2C2 is an external dependency and is not vendored by this repository.

    if nargin < 2 || strlength(string(solverName)) == 0
        solverName = "glpk";
    end

    model = flutor.validateModel(model);

    if exist("F2C2", "file") ~= 2
        error( ...
            "flutor:F2C2DependencyMissing", ...
            ["F2C2 was not found on the MATLAB path. Install the F2C2 " ...
             "implementation used for the analysis, or run with a precomputed " ...
             "full-coupling matrix."]);
    end

    f2c2Model = struct( ...
        "stoichiometricMatrix", model.S, ...
        "Reactions", {model.rxns}, ...
        "Metabolites", {model.mets}, ...
        "reversibilityVector", double(model.lb < 0), ...
        "lb", model.lb, ...
        "ub", model.ub);

    [rawTable, ~] = F2C2(char(string(solverName)), f2c2Model);
    rawTable = double(rawTable);

    nReactions = numel(model.rxns);
    if ~isequal(size(rawTable), [nReactions, nReactions])
        error( ...
            "flutor:InvalidF2C2Output", ...
            "F2C2 returned a matrix with unexpected dimensions.");
    end
    if any(~isfinite(rawTable), "all")
        error("flutor:InvalidF2C2Output", "F2C2 returned non-finite values.");
    end

    fullCoupled = rawTable == 1;
    fullCoupled = fullCoupled | fullCoupled.';
    fullCoupled(1:nReactions + 1:end) = true;

    diagnostics = struct( ...
        "solver", string(solverName), ...
        "nFullCoupledPairs", nnz(triu(fullCoupled, 1)));
end
