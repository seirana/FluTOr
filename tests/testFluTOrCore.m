function tests = testFluTOrCore
%TESTFLUTORCORE Solver-independent tests for maintained FluTOr utilities.

    tests = functiontests(localfunctions);
end

function testValidateModelAddsIndices(testCase)
    model = flutor.validateModel(localToyModel());
    verifyEqual(testCase, model.rxnNumber, (1:3)');
    verifyEqual(testCase, model.metNumber, (1:2)');
end

function testApplyReactionBoundsByID(testCase)
    model = localToyModel();
    changes = table( ...
        ["R1"; "R3"], ...
        [-5; 2], ...
        [-1; 2], ...
        'VariableNames', ...
        {'ReactionID', 'LowerBound', 'UpperBound'});

    changed = flutor.applyReactionBounds(model, changes);

    verifyEqual(testCase, changed.lb, [-5; 0; 2]);
    verifyEqual(testCase, changed.ub, [-1; 10; 2]);
end

function testCombineCouplingEvidenceUsesIntersection(testCase)
    qfca = logical([ ...
        1, 1, 0, 0; ...
        1, 1, 0, 0; ...
        0, 0, 1, 1; ...
        0, 0, 1, 1]);
    f2c2 = logical([ ...
        1, 1, 0, 0; ...
        1, 1, 0, 0; ...
        0, 0, 1, 0; ...
        0, 0, 0, 1]);

    coupling = flutor.combineFullCouplingEvidence(qfca, f2c2);

    verifyTrue(testCase, coupling.consensus(1, 2));
    verifyFalse(testCase, coupling.consensus(3, 4));
    verifyEqual(testCase, coupling.nConsensusPairs, 1);
end

function testExpandCoupledAlternativesIncludesSelf(testCase)
    matrix = eye(4) > 0;
    matrix(1, 3) = true;
    matrix(3, 1) = true;
    matrix(2, 4) = true;
    matrix(4, 2) = true;

    combinations = flutor.expandCoupledAlternatives([1, 2], matrix);
    combinations = sortrows(combinations);

    expected = sortrows([ ...
        1, 2; ...
        1, 4; ...
        2, 3; ...
        3, 4]);

    verifyEqual(testCase, combinations, expected);
end

function testQfcaTrivialFullCoupling(testCase)
    model = struct();
    model.S = [1, -1];
    model.rxns = {"R1"; "R2"};
    model.mets = {"M1"};
    model.lb = [0; 0];
    model.ub = [10; 10];
    model.c = [0; 0];

    [matrix, diagnostics] = ...
        flutor.computeQfcaFullCoupling(model, 1e-10);

    verifyTrue(testCase, matrix(1, 2));
    verifyTrue(testCase, matrix(2, 1));
    verifyEqual(testCase, diagnostics.nOriginalReactions, 2);
    verifyGreaterThanOrEqual(testCase, diagnostics.nTrivialMerges, 1);
end

function testUnknownReactionBoundFails(testCase)
    model = localToyModel();
    changes = table( ...
        "missing", ...
        0, ...
        1, ...
        'VariableNames', ...
        {'ReactionID', 'LowerBound', 'UpperBound'});

    verifyError( ...
        testCase, ...
        @() flutor.applyReactionBounds(model, changes), ...
        "flutor:UnknownReaction");
end

function model = localToyModel()
    model = struct();
    model.S = [ ...
        1, -1, 0; ...
        0, 1, -1];
    model.rxns = {"R1"; "R2"; "R3"};
    model.mets = {"M1"; "M2"};
    model.lb = [0; 0; 0];
    model.ub = [10; 10; 10];
    model.c = zeros(3, 1);
end
