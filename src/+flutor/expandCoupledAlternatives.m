function combinations = expandCoupledAlternatives(selectedIndices, fullCoupledMatrix)
%EXPANDCOUPLEDALTERNATIVES Expand selected reactions through coupling groups.

    selected = double(selectedIndices(:)');
    matrix = logical(fullCoupledMatrix);

    if size(matrix, 1) ~= size(matrix, 2)
        error("flutor:InvalidCouplingMatrix", "fullCoupledMatrix must be square.");
    end
    n = size(matrix, 1);
    if any(selected < 1) || any(selected > n) || any(mod(selected, 1) ~= 0)
        error("flutor:InvalidReactionIndex", "selectedIndices contains invalid indices.");
    end
    if isempty(selected)
        combinations = zeros(0, 0);
        return;
    end

    choices = cell(1, numel(selected));
    for i = 1:numel(selected)
        choices{i} = unique( ...
            [selected(i), find(matrix(selected(i), :))], ...
            "stable");
    end

    grids = cell(1, numel(choices));
    [grids{:}] = ndgrid(choices{:});
    combinations = zeros(numel(grids{1}), numel(choices));

    for i = 1:numel(grids)
        combinations(:, i) = grids{i}(:);
    end

    valid = false(size(combinations, 1), 1);
    for row = 1:size(combinations, 1)
        combinations(row, :) = sort(combinations(row, :));
        valid(row) = numel(unique(combinations(row, :))) == numel(selected);
    end
    combinations = unique(combinations(valid, :), "rows", "stable");
end
