function [model, diagnostics] = computeFluxRanges(model, biomassReaction, options)
%COMPUTEFLUXRANGES Compute FBA/FVA with checked solver status.
%
%   The function uses Optimization Toolbox linprog so additional linear
%   inequalities can be represented explicitly. It never converts an
%   infeasible optimization into a zero flux range.
%
%   OPTIONS fields:
%     biomassFraction    fraction of optimum biomass to enforce, [] disables
%     Aineq, bineq       additional linear inequalities
%     tolerance          zero/block tolerance
%     convertIrreversible convert with COBRA Toolbox after FVA (default true)

    if nargin < 3 || isempty(options)
        options = struct();
    end
    options = localOptions(options);

    if exist("linprog", "file") ~= 2
        error( ...
            "flutor:OptimizationToolboxMissing", ...
            "Optimization Toolbox function linprog is required.");
    end

    model = flutor.validateModel(model);
    biomassIndex = localReactionIndex(model, biomassReaction);

    nReactions = numel(model.rxns);
    objective = zeros(nReactions, 1);
    objective(biomassIndex) = 1;
    model.c = objective;

    Aeq = model.S;
    beq = zeros(size(Aeq, 1), 1);
    lower = model.lb;
    upper = model.ub;

    if isempty(options.Aineq)
        Aineq = zeros(0, nReactions);
        bineq = zeros(0, 1);
    else
        Aineq = double(options.Aineq);
        bineq = double(options.bineq(:));
        if size(Aineq, 2) ~= nReactions || size(Aineq, 1) ~= numel(bineq)
            error( ...
                "flutor:InvalidInequalityDimensions", ...
                "Aineq/bineq dimensions are incompatible with the model.");
        end
    end

    linprogOptions = optimoptions("linprog", "Display", "none");

    optimumBiomass = NaN;
    if ~isempty(options.biomassFraction)
        maximizeObjective = -objective;
        [~, objectiveValue, exitFlag, solverOutput] = linprog( ...
            maximizeObjective, ...
            Aineq, ...
            bineq, ...
            Aeq, ...
            beq, ...
            lower, ...
            upper, ...
            linprogOptions);
        localRequireOptimal(exitFlag, solverOutput, "biomass maximization");

        optimumBiomass = -objectiveValue;
        if ~isfinite(optimumBiomass)
            error("flutor:OptimizationFailed", "Biomass optimum is non-finite.");
        end

        lower(biomassIndex) = options.biomassFraction * optimumBiomass;
        if lower(biomassIndex) > upper(biomassIndex)
            error( ...
                "flutor:InfeasibleBiomassFraction", ...
                "Requested biomass fraction exceeds the biomass upper bound.");
        end
    end

    minima = zeros(nReactions, 1);
    maxima = zeros(nReactions, 1);

    for i = 1:nReactions
        unitObjective = zeros(nReactions, 1);
        unitObjective(i) = 1;

        [~, minimumValue, exitFlagMin, outputMin] = linprog( ...
            unitObjective, Aineq, bineq, Aeq, beq, lower, upper, linprogOptions);
        localRequireOptimal(exitFlagMin, outputMin, sprintf("minimum of reaction %s", model.rxns{i}));

        [~, negativeMaximum, exitFlagMax, outputMax] = linprog( ...
            -unitObjective, Aineq, bineq, Aeq, beq, lower, upper, linprogOptions);
        localRequireOptimal(exitFlagMax, outputMax, sprintf("maximum of reaction %s", model.rxns{i}));

        minima(i) = minimumValue;
        maxima(i) = -negativeMaximum;

        if minima(i) > maxima(i) + options.tolerance
            error( ...
                "flutor:InconsistentFluxRange", ...
                "Minimum exceeds maximum for reaction %s.", ...
                model.rxns{i});
        end
    end

    model.lb = minima;
    model.ub = maxima;

    nBefore = nReactions;
    if options.convertIrreversible
        if exist("convertToIrreversible", "file") ~= 2
            error( ...
                "flutor:CobraDependencyMissing", ...
                "COBRA Toolbox function convertToIrreversible was not found.");
        end
        model = convertToIrreversible(model);
        model = flutor.validateModel(model);
    end

    model = localCanonicalizeNegativeIrreversible(model, options.tolerance);
    [model, removed] = localRemoveBlocked(model, options.tolerance);

    diagnostics = struct( ...
        "biomassReaction", string(model.rxns(localSafeBiomassIndex(model, biomassReaction))), ...
        "biomassFraction", options.biomassFraction, ...
        "optimumBiomass", optimumBiomass, ...
        "nReactionsBeforeIrreversibleConversion", nBefore, ...
        "nReactionsAfterPreprocessing", numel(model.rxns), ...
        "nBlockedReactionsRemoved", removed.nReactions, ...
        "nDeadEndMetabolitesRemoved", removed.nMetabolites, ...
        "tolerance", options.tolerance);
end

function options = localOptions(options)
    defaults = struct( ...
        "biomassFraction", [], ...
        "Aineq", [], ...
        "bineq", [], ...
        "tolerance", 1e-5, ...
        "convertIrreversible", true);

    names = fieldnames(defaults);
    for i = 1:numel(names)
        name = names{i};
        if ~isfield(options, name)
            options.(name) = defaults.(name);
        end
    end

    if ~isempty(options.biomassFraction)
        validateattributes(options.biomassFraction, {"numeric"}, {"scalar", "real", "positive", "<=", 1});
    end
    validateattributes(options.tolerance, {"numeric"}, {"scalar", "real", "positive", "finite"});
    options.convertIrreversible = logical(options.convertIrreversible);
end

function index = localReactionIndex(model, reaction)
    if isnumeric(reaction)
        if ~isscalar(reaction) || reaction < 1 || reaction > numel(model.rxns) || mod(reaction, 1) ~= 0
            error("flutor:UnknownBiomassReaction", "Invalid biomass reaction index.");
        end
        index = double(reaction);
        return;
    end

    matches = find(string(model.rxns) == string(reaction));
    if numel(matches) ~= 1
        error( ...
            "flutor:UnknownBiomassReaction", ...
            "Biomass reaction %s must identify exactly one reaction.", ...
            string(reaction));
    end
    index = matches;
end

function localRequireOptimal(exitFlag, output, context)
    if exitFlag ~= 1
        message = "";
        if isstruct(output) && isfield(output, "message")
            message = string(output.message);
        end
        error( ...
            "flutor:OptimizationFailed", ...
            "linprog failed during %s (exitFlag=%d). %s", ...
            context, ...
            exitFlag, ...
            message);
    end
end

function model = localCanonicalizeNegativeIrreversible(model, tolerance)
    flip = model.lb < -tolerance & model.ub <= tolerance;
    oldLower = model.lb(flip);
    oldUpper = model.ub(flip);
    model.lb(flip) = -oldUpper;
    model.ub(flip) = -oldLower;
    model.S(:, flip) = -model.S(:, flip);
end

function [model, removed] = localRemoveBlocked(model, tolerance)
    zeroColumns = full(sum(abs(model.S), 1))' <= tolerance;
    zeroRanges = abs(model.lb) <= tolerance & abs(model.ub) <= tolerance;
    reactionMask = zeroColumns | zeroRanges;

    nReactions = nnz(reactionMask);
    nMetabolites = 0;

    if any(reactionMask)
        if exist("removeRxns", "file") ~= 2
            error( ...
                "flutor:CobraDependencyMissing", ...
                "COBRA Toolbox function removeRxns is required.");
        end

        originalRxns = string(model.rxns);
        originalNumbers = model.rxnNumber;
        model = removeRxns(model, model.rxns(reactionMask));
        model = flutor.validateModel(model);
        [found, indices] = ismember(string(model.rxns), originalRxns);
        if ~all(found)
            error("flutor:ReactionIdentityLost", "Could not restore reaction-number mapping.");
        end
        model.rxnNumber = originalNumbers(indices);
    end

    deadMetabolites = full(sum(abs(model.S), 2)) <= tolerance;
    if any(deadMetabolites) && exist("removeMetabolites", "file") == 2
        nMetabolites = nnz(deadMetabolites);
        model = removeMetabolites(model, model.mets(deadMetabolites), false);
        model = flutor.validateModel(model);
    end

    removed = struct("nReactions", nReactions, "nMetabolites", nMetabolites);
end

function index = localSafeBiomassIndex(model, biomassReaction)
    if isnumeric(biomassReaction)
        if isfield(model, "rxnNumber")
            match = find(model.rxnNumber == biomassReaction, 1);
            if ~isempty(match)
                index = match;
                return;
            end
        end
        index = 1;
        return;
    end

    match = find(string(model.rxns) == string(biomassReaction), 1);
    if isempty(match)
        index = 1;
    else
        index = match;
    end
end
