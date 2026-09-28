function coupling = combineFullCouplingEvidence(qfcaMatrix, f2c2Matrix)
%COMBINEFULLCOUPLINGEVIDENCE Intersect full-coupling calls from two methods.

    qfca = logical(qfcaMatrix);
    f2c2 = logical(f2c2Matrix);

    if ~isequal(size(qfca), size(f2c2)) ...
            || size(qfca, 1) ~= size(qfca, 2)
        error( ...
            "flutor:CouplingDimensionMismatch", ...
            "QFCA and F2C2 coupling matrices must be square and equal-sized.");
    end

    n = size(qfca, 1);
    qfca = qfca | qfca.';
    f2c2 = f2c2 | f2c2.';
    qfca(1:n + 1:end) = true;
    f2c2(1:n + 1:end) = true;

    consensus = qfca & f2c2;
    consensus = localTransitiveClosure(consensus);
    consensus(1:n + 1:end) = true;

    coupling = struct( ...
        "qfca", qfca, ...
        "f2c2", f2c2, ...
        "consensus", consensus, ...
        "groups", {localGroups(consensus)}, ...
        "nConsensusPairs", nnz(triu(consensus, 1)));
end

function closure = localTransitiveClosure(matrix)
    closure = logical(matrix);
    n = size(closure, 1);

    for k = 1:n
        closure = closure | (closure(:, k) & closure(k, :));
    end
    closure = closure | closure.';
end

function groups = localGroups(matrix)
    n = size(matrix, 1);
    assigned = false(n, 1);
    groups = {};

    for i = 1:n
        if assigned(i)
            continue;
        end

        members = find(matrix(i, :));
        members = unique([i, members], "stable");
        assigned(members) = true;
        groups{end + 1, 1} = members(:); %#ok<AGROW>
    end
end
