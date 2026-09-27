function [fullCoupled, diagnostics] = computeQfcaFullCoupling(model, tolerance)
%COMPUTEQFCAFULLCOUPLING Compute QFCA fully-coupled equivalence classes.
%
%   This maintained implementation isolates the fully-coupled portion of the
%   modified QFCA kernel used by FluTOr. Directional-coupling codes are not
%   computed because FluTOr's downstream consensus step uses only full
%   coupling evidence.
%
%   Two defects from the historical implementation are avoided:
%     * the reduced reaction count is updated only when a merge occurs;
%     * no uninitialized LP-counter state is referenced.

    if nargin < 2 || isempty(tolerance)
        tolerance = 1e-10;
    end
    validateattributes(tolerance, {'numeric'}, {'scalar', 'real', 'positive', 'finite'});

    model = flutor.validateModel(model);
    nOriginal = numel(model.rxns);

    reduced = struct();
    reduced.S = sparse(model.S);
    reduced.rev = model.lb < -tolerance;
    reduced.lb = model.lb;
    reduced.ub = model.ub;

    % Remove duplicate metabolite constraints, matching publication-kernel behavior.
    [reduced.S, ~, ~] = unique(reduced.S, "rows", "stable");

    % Aggregate identical reaction columns with identical reversibility.
    signature = [reduced.S.', sparse(double(reduced.rev))];
    [~, representativeIndices, duplicateMap] = unique( ...
        signature, ...
        "rows", ...
        "stable");

    reduced.S = reduced.S(:, representativeIndices);
    reduced.rev = reduced.rev(representativeIndices);
    reduced.lb = reduced.lb(representativeIndices);
    reduced.ub = reduced.ub(representativeIndices);

    representativeReaction = representativeIndices(:);
    fullCouplingRepresentative = representativeReaction(duplicateMap);

    % Trivial full couplings induced by two-reaction metabolite balances.
    merged = true;
    trivialMergeCount = 0;

    while merged
        merged = false;
        [nMetabolites, ~] = size(reduced.S);

        for metaboliteIndex = nMetabolites:-1:1
            nonzeroColumns = find(reduced.S(metaboliteIndex, :));
            if numel(nonzeroColumns) ~= 2
                continue;
            end

            i = nonzeroColumns(1);
            j = nonzeroColumns(2);

            fixedI = abs(reduced.ub(i) - reduced.lb(i)) <= tolerance;
            fixedJ = abs(reduced.ub(j) - reduced.lb(j)) <= tolerance;
            if fixedI ~= fixedJ
                continue;
            end

            denominator = reduced.S(metaboliteIndex, j);
            if abs(denominator) <= tolerance
                continue;
            end

            coefficient = -reduced.S(metaboliteIndex, i) / denominator;
            reduced = localMerge(reduced, i, j, coefficient);

            oldRepresentative = representativeReaction(j);
            newRepresentative = representativeReaction(i);
            fullCouplingRepresentative( ...
                fullCouplingRepresentative == oldRepresentative) = newRepresentative;
            representativeReaction(j) = [];

            trivialMergeCount = trivialMergeCount + 1;
            merged = true;
            break;
        end
    end

    % Remaining full couplings from normalized null-space collinearity.
    qrMergeCount = 0;
    if size(reduced.S, 2) >= 2
        [Q, R] = qr(reduced.S.');
        numericalRankTolerance = max( ...
            tolerance, ...
            norm(reduced.S, "fro") * eps(class(full(reduced.S))));
        diagonal = abs(diag(R));
        rankS = sum(diagonal > numericalRankTolerance);
        nReduced = size(reduced.S, 2);

        if rankS < nReduced
            Z = Q(:, rankS + 1:nReduced);
            gram = Z * Z.';
            diagonalGram = diag(gram);
            scale = zeros(size(diagonalGram));
            positive = diagonalGram > numericalRankTolerance^2;
            scale(positive) = 1 ./ sqrt(diagonalGram(positive));
            normalized = (scale .* gram) .* scale.';

            for i = nReduced:-1:2
                candidates = abs(normalized(i, 1:i-1));
                [maximumCorrelation, j] = max(candidates);

                if isempty(maximumCorrelation) ...
                        || maximumCorrelation <= 1 - numericalRankTolerance
                    continue;
                end

                fixedI = abs(reduced.ub(i) - reduced.lb(i)) <= tolerance;
                fixedJ = abs(reduced.ub(j) - reduced.lb(j)) <= tolerance;
                if fixedI ~= fixedJ
                    continue;
                end

                if scale(i) == 0 || scale(j) == 0
                    continue;
                end

                coefficient = sign(normalized(i, j)) * scale(j) / scale(i);
                reduced = localMerge(reduced, j, i, coefficient);

                oldRepresentative = representativeReaction(i);
                newRepresentative = representativeReaction(j);
                fullCouplingRepresentative( ...
                    fullCouplingRepresentative == oldRepresentative) = newRepresentative;
                representativeReaction(i) = [];

                normalized(i, :) = [];
                normalized(:, i) = [];
                scale(i) = [];
                qrMergeCount = qrMergeCount + 1;
            end
        end
    end

    fullCoupled = ...
        fullCouplingRepresentative(:) == fullCouplingRepresentative(:).';
    fullCoupled = logical(fullCoupled);
    fullCoupled(1:nOriginal + 1:end) = true;

    diagnostics = struct( ...
        "nOriginalReactions", nOriginal, ...
        "nInitialIsozymeGroups", numel(representativeIndices), ...
        "nTrivialMerges", trivialMergeCount, ...
        "nQrMerges", qrMergeCount, ...
        "nFinalFullCouplingGroups", numel(unique(fullCouplingRepresentative)), ...
        "tolerance", tolerance);
end

function reduced = localMerge(reduced, keepIndex, removeIndex, coefficient)
    reduced.S(:, keepIndex) = ...
        reduced.S(:, keepIndex) + coefficient * reduced.S(:, removeIndex);
    reduced.S(:, removeIndex) = [];
    reduced.lb(removeIndex) = [];
    reduced.ub(removeIndex) = [];
    reduced.rev(removeIndex) = [];
end
