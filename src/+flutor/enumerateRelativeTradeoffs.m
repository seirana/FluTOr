function result = enumerateRelativeTradeoffs(model, biomassReaction, coupling, options)
%ENUMERATERELATIVETRADEOFFS Enumerate relative trade-offs with checked MILP.
%
%   The maintained formulation preserves the publication-era decision
%   structure while returning structured results instead of saving MAT files
%   during each solver iteration.

    if nargin < 4 || isempty(options)
        options = struct();
    end
    options = localOptions(options);

    if exist("intlinprog", "file") ~= 2
        error( ...
            "flutor:OptimizationToolboxMissing", ...
            "Optimization Toolbox function intlinprog is required.");
    end

    model = flutor.validateModel(model);
    nReactions = numel(model.rxns);

    if ~isstruct(coupling) || ~isfield(coupling, "consensus")
        error("flutor:InvalidCoupling", "coupling.consensus is required.");
    end
    fullCoupled = logical(coupling.consensus);
    if ~isequal(size(fullCoupled), [nReactions, nReactions])
        error("flutor:InvalidCoupling", "coupling.consensus has incompatible dimensions.");
    end

    biomassIndex = localReactionIndex(model, biomassReaction);
    biomassGroup = find(fullCoupled(biomassIndex, :));
    biomassGroup = unique([biomassIndex, biomassGroup]);

    fixedMask = abs(model.ub - model.lb) <= options.tolerance;
    fixedIndices = setdiff(find(fixedMask), biomassGroup);

    variableMask = ~fixedMask;
    variableIndices = setdiff(find(variableMask), biomassGroup);
    variableIndices = sort(variableIndices(:));

    if isempty(variableIndices)
        result = localEmptyResult(model, biomassGroup, variableIndices, options);
        return;
    end

    nMetabolites = size(model.S, 1);
    nVariables = numel(variableIndices);

    % Decision vector: [k(m); alpha(r); u(v); s(v)]
    alphaStart = nMetabolites + 1;
    uStart = nMetabolites + nReactions + 1;
    sStart = uStart + nVariables;
    totalVariables = nMetabolites + nReactions + 2 * nVariables;

    Aeq = [ ...
        sparse(model.S'), ...
        -speye(nReactions), ...
        sparse(nReactions, nVariables), ...
        sparse(nReactions, nVariables)];
    beq = zeros(nReactions, 1);

    Ain = sparse(0, totalVariables);
    bin = zeros(0, 1);

    % At least two selected variable reactions.
    row = sparse(1, totalVariables);
    row(sStart:(sStart + nVariables - 1)) = -1;
    Ain(end + 1, :) = row;
    bin(end + 1, 1) = -2;

    % Positive aggregate alpha contribution from the biomass coupling group.
    row = sparse(1, totalVariables);
    row(alphaStart + biomassGroup - 1) = -1;
    Ain(end + 1, :) = row;
    bin(end + 1, 1) = -1;

    selector = sparse(nVariables, nReactions);
    for i = 1:nVariables
        selector(i, variableIndices(i)) = 1;
    end

    Ain = [ ...
        Ain; ...
        sparse(nVariables, nMetabolites), selector, -speye(nVariables), sparse(nVariables, nVariables); ...
        sparse(nVariables, nMetabolites), -selector, -speye(nVariables), sparse(nVariables, nVariables); ...
        sparse(nVariables, nMetabolites), selector, sparse(nVariables, nVariables), options.bigM * speye(nVariables); ...
        sparse(nVariables, nMetabolites), -selector, sparse(nVariables, nVariables), -options.bigM * speye(nVariables)];
    bin = [ ...
        bin; ...
        zeros(nVariables, 1); ...
        zeros(nVariables, 1); ...
        (options.bigM - 1) * ones(nVariables, 1); ...
        zeros(nVariables, 1)];

    % Prevent simultaneous selection of fully-coupled alternatives outside
    % the biomass group.
    groups = localCouplingGroups(fullCoupled);
    for groupIndex = 1:numel(groups)
        members = setdiff(groups{groupIndex}, biomassGroup);
        localMembers = find(ismember(variableIndices, members));
        if numel(localMembers) <= 1
            continue;
        end
        row = sparse(1, totalVariables);
        row(sStart + localMembers - 1) = 1;
        Ain(end + 1, :) = row;
        bin(end + 1, 1) = 1;
    end

    alphaLower = zeros(nReactions, 1);
    alphaUpper = zeros(nReactions, 1);

    alphaLower(biomassGroup) = 0;
    alphaUpper(biomassGroup) = Inf;

    alphaLower(fixedIndices) = -Inf;
    alphaUpper(fixedIndices) = Inf;

    alphaLower(variableIndices) = -options.coefficientLimit;
    alphaUpper(variableIndices) = 0;

    lb = [ ...
        -Inf(nMetabolites, 1); ...
        alphaLower; ...
        zeros(nVariables, 1); ...
        zeros(nVariables, 1)];
    ub = [ ...
        Inf(nMetabolites, 1); ...
        alphaUpper; ...
        options.coefficientLimit * ones(nVariables, 1); ...
        ones(nVariables, 1)];

    objective = zeros(totalVariables, 1);
    objective(sStart:(sStart + nVariables - 1)) = 1;

    integerVariables = sStart:(sStart + nVariables - 1);
    solverOptions = optimoptions( ...
        "intlinprog", ...
        "Display", options.display, ...
        "IntegerPreprocess", "none", ...
        "LPPreprocess", "none", ...
        "MaxTime", options.maxTimeSeconds);

    supportRows = false(0, nReactions);
    solverRuns = 0;

    while size(supportRows, 1) < options.maxSolutions
        [solution, objectiveValue, exitFlag, solverOutput] = intlinprog( ...
            objective, ...
            integerVariables, ...
            Ain, ...
            bin, ...
            Aeq, ...
            beq, ...
            lb, ...
            ub, ...
            solverOptions);
        solverRuns = solverRuns + 1;

        if exitFlag == -2
            break;
        end
        if exitFlag <= 0 || isempty(solution)
            error( ...
                "flutor:TradeoffSolverFailed", ...
                "intlinprog failed (exitFlag=%d): %s", ...
                exitFlag, ...
                localSolverMessage(solverOutput));
        end

        selectedBinary = solution(sStart:(sStart + nVariables - 1)) >= 0.5;
        selected = variableIndices(selectedBinary);
        degree = numel(selected);

        if degree < 2
            error( ...
                "flutor:UnexpectedMILPSupport", ...
                "Solver returned a trade-off support with degree %d.", ...
                degree);
        end

        alternatives = flutor.expandCoupledAlternatives( ...
            selected, ...
            fullCoupled);
        if isempty(alternatives)
            alternatives = selected(:)';
        end

        for alternativeIndex = 1:size(alternatives, 1)
            support = false(1, nReactions);
            support(alternatives(alternativeIndex, :)) = true;

            if ~any(all(supportRows == support, 2))
                supportRows(end + 1, :) = support; %#ok<AGROW>
            end

            localSelected = find(ismember(variableIndices, alternatives(alternativeIndex, :)));
            if ~isempty(localSelected)
                row = sparse(1, totalVariables);
                row(sStart + localSelected - 1) = 1;
                Ain(end + 1, :) = row;
                bin(end + 1, 1) = numel(localSelected) - 1;
            end
        end

        % Exclude the exact binary support regardless of expansion overlap.
        exactLocal = find(selectedBinary);
        row = sparse(1, totalVariables);
        row(sStart + exactLocal - 1) = 1;
        Ain(end + 1, :) = row;
        bin(end + 1, 1) = numel(exactLocal) - 1;

        if ~isfinite(objectiveValue)
            error("flutor:TradeoffSolverFailed", "MILP objective was non-finite.");
        end
    end

    result = struct( ...
        "tradeoffs", localTradeoffTable(model, biomassGroup, supportRows), ...
        "supportMatrix", supportRows, ...
        "biomassGroup", biomassGroup(:), ...
        "variableReactionIndices", variableIndices, ...
        "solverRuns", solverRuns, ...
        "options", options);
end

function options = localOptions(options)
    defaults = struct( ...
        "tolerance", 1e-5, ...
        "coefficientLimit", 1000, ...
        "bigM", 1001, ...
        "maxTimeSeconds", 3600, ...
        "maxSolutions", Inf, ...
        "display", "off");

    names = fieldnames(defaults);
    for i = 1:numel(names)
        name = names{i};
        if ~isfield(options, name)
            options.(name) = defaults.(name);
        end
    end

    validateattributes(options.tolerance, {"numeric"}, {"scalar", "real", "positive", "finite"});
    validateattributes(options.coefficientLimit, {"numeric"}, {"scalar", "real", "positive", "finite"});
    validateattributes(options.bigM, {"numeric"}, {"scalar", "real", "positive", "finite"});
    validateattributes(options.maxTimeSeconds, {"numeric"}, {"scalar", "real", "positive", "finite"});
    validateattributes(options.maxSolutions, {"numeric"}, {"scalar", "real", "positive"});
    options.display = char(string(options.display));
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

function groups = localCouplingGroups(matrix)
    n = size(matrix, 1);
    assigned = false(n, 1);
    groups = {};

    for i = 1:n
        if assigned(i)
            continue;
        end
        members = unique([i, find(matrix(i, :))], "stable");
        assigned(members) = true;
        groups{end + 1, 1} = members(:); %#ok<AGROW>
    end
end

function message = localSolverMessage(output)
    message = "No solver message available.";
    if isstruct(output) && isfield(output, "message") && ~isempty(output.message)
        message = string(output.message);
    end
end

function tableOut = localTradeoffTable(model, biomassGroup, supportRows)
    nTradeoffs = size(supportRows, 1);
    tradeoffID = (1:nTradeoffs)';
    degree = zeros(nTradeoffs, 1);
    reactionNumbers = strings(nTradeoffs, 1);
    reactionIDs = strings(nTradeoffs, 1);

    for i = 1:nTradeoffs
        members = find(supportRows(i, :));
        degree(i) = numel(members);
        reactionNumbers(i) = strjoin(string(model.rxnNumber(members)), ";");
        reactionIDs(i) = strjoin(string(model.rxns(members)), ";");
    end

    biomassIDs = strjoin(string(model.rxns(biomassGroup)), ";");
    biomassColumn = repmat(biomassIDs, nTradeoffs, 1);

    tableOut = table( ...
        tradeoffID, ...
        degree, ...
        reactionNumbers, ...
        reactionIDs, ...
        biomassColumn, ...
        "VariableNames", ...
        {"TradeoffID", "Degree", "ReactionNumbers", "ReactionIDs", "BiomassCouplingGroup"});
end

function result = localEmptyResult(model, biomassGroup, variableIndices, options)
    result = struct( ...
        "tradeoffs", table( ...
            zeros(0, 1), zeros(0, 1), strings(0, 1), strings(0, 1), strings(0, 1), ...
            "VariableNames", ...
            {"TradeoffID", "Degree", "ReactionNumbers", "ReactionIDs", "BiomassCouplingGroup"}), ...
        "supportMatrix", false(0, numel(model.rxns)), ...
        "biomassGroup", biomassGroup(:), ...
        "variableReactionIndices", variableIndices, ...
        "solverRuns", 0, ...
        "options", options);
end
