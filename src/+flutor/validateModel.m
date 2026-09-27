function model = validateModel(model)
%VALIDATEMODEL Validate and normalize the model fields required by FluTOr.

    if ~isstruct(model)
        error("flutor:InvalidModel", "Model must be a struct.");
    end

    required = ["S", "rxns", "mets", "lb", "ub", "c"];
    missing = required(~isfield(model, cellstr(required)));
    if ~isempty(missing)
        error( ...
            "flutor:MissingModelFields", ...
            "Model is missing required fields: %s", ...
            strjoin(missing, ", "));
    end

    if (~isnumeric(model.S) && ~issparse(model.S)) || ndims(model.S) ~= 2
        error("flutor:InvalidStoichiometry", "model.S must be a two-dimensional numeric matrix.");
    end
    if any(~isfinite(nonzeros(model.S)))
        error("flutor:InvalidStoichiometry", "model.S contains non-finite values.");
    end

    nMetabolites = size(model.S, 1);
    nReactions = size(model.S, 2);

    if numel(model.rxns) ~= nReactions || numel(model.mets) ~= nMetabolites
        error( ...
            "flutor:DimensionMismatch", ...
            "Reaction/metabolite identifiers must match model.S dimensions.");
    end

    model.rxns = cellstr(string(model.rxns(:)));
    model.mets = cellstr(string(model.mets(:)));
    model.lb = double(model.lb(:));
    model.ub = double(model.ub(:));
    model.c = double(model.c(:));

    if numel(model.lb) ~= nReactions ...
            || numel(model.ub) ~= nReactions ...
            || numel(model.c) ~= nReactions
        error( ...
            "flutor:DimensionMismatch", ...
            "lb, ub, and c must contain one value per reaction.");
    end

    if any(~isfinite(model.lb)) || any(~isfinite(model.ub)) || any(~isfinite(model.c))
        error("flutor:NonFiniteModelValues", "Bounds/objective contain non-finite values.");
    end
    if any(model.lb > model.ub)
        error("flutor:InvalidBounds", "Every lower bound must be <= its upper bound.");
    end

    if ~isfield(model, "rxnNumber") || numel(model.rxnNumber) ~= nReactions
        model.rxnNumber = (1:nReactions)';
    else
        model.rxnNumber = double(model.rxnNumber(:));
    end

    if ~isfield(model, "metNumber") || numel(model.metNumber) ~= nMetabolites
        model.metNumber = (1:nMetabolites)';
    else
        model.metNumber = double(model.metNumber(:));
    end
end
