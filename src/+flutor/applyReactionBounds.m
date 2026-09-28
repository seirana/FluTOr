function model = applyReactionBounds(model, changes)
%APPLYREACTIONBOUNDS Apply an explicit table of reaction-bound changes.
%
%   CHANGES must be a table with variables:
%     ReactionID, LowerBound, UpperBound

    model = flutor.validateModel(model);

    if isempty(changes)
        return;
    end
    if ~istable(changes)
        error("flutor:InvalidBoundChanges", "changes must be a table.");
    end

    required = ["ReactionID", "LowerBound", "UpperBound"];
    missing = required(~ismember(required, string(changes.Properties.VariableNames)));
    if ~isempty(missing)
        error( ...
            "flutor:InvalidBoundChanges", ...
            "Missing table variables: %s", ...
            strjoin(missing, ", "));
    end

    reactionIDs = string(changes.ReactionID);
    lower = double(changes.LowerBound);
    upper = double(changes.UpperBound);

    if any(~isfinite(lower)) || any(~isfinite(upper))
        error("flutor:InvalidBoundChanges", "Bounds must be finite.");
    end
    if any(lower > upper)
        error("flutor:InvalidBoundChanges", "Each lower bound must be <= its upper bound.");
    end
    if numel(unique(reactionIDs)) ~= numel(reactionIDs)
        error("flutor:DuplicateBoundChange", "ReactionID values must be unique.");
    end

    [found, indices] = ismember(reactionIDs, string(model.rxns));
    if any(~found)
        error( ...
            "flutor:UnknownReaction", ...
            "Unknown reaction IDs: %s", ...
            strjoin(reactionIDs(~found), ", "));
    end

    model.lb(indices) = lower;
    model.ub(indices) = upper;
end
