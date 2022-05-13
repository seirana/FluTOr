function QFCA = mergeFullyCoupled(QFCA, i, j, c)
% edited by Seirana
%
% mergeFullyCoupled merges the fully coupled pair of reactions (i, j)
%
% USAGE:
%    [S, rev, rxns, lb, ub, rxnN, FC] = mergeFullyCoupled(S, rev, rxns, i, j, c, lb, ub, rxnN, FC)
%
%     INPUTS:
%         S:            the associated sparse stoichiometric matrix
%         rev:          the 0-1 vector with 1's corresponding to the reversible reactions
%         rxns:         the cell array of reaction abbreviations
%         i:            the reaction which the other one will merge with and is not removed
%         j:            the reaction which will be merged into the other reaction
%         c:            the full coupling coefficient
%         lb:           the doulbe array of reaction flux lower bound
%         ub:           the doulbe array of reaction flux upper bound
%         rxnN:         the rxn numbers in original model
%         FC: the FC array
%
%     OUTPUTS:
%         S:            the reduced sparse stoichiometric matrix
%         rev:          the reduced reversibility vector
%         rxns:         the reduced reaction abbreviations
%         lb:           the reduced lower bound vector
%         ub:           the reduced upper bound vector
%         rxnN:         the rxn numbers in original model
%         FC: the FC array
%
% .. Authors:
%        - Mojtaba Tefagh, Stephen P. Boyd, 2019, Stanford University

QFCA.S(:, i) = QFCA.S(:, i) + c*QFCA.S(:, j);
QFCA.S(:, j) = [];
QFCA.lb(j) = [];
QFCA.ub(j) = [];
QFCA.rev(j) = [];
QFCA.rxns(i) = {strjoin([QFCA.rxns(i), QFCA.rxns(j)], ', ')};
QFCA.rxns(j) = [];
QFCA.rxnN(i) = {strjoin([QFCA.rxnN(i), QFCA.rxnN(j)], ', ')};
QFCA.rxnN(j) = [];
QFCA.FC(i) = {strjoin([QFCA.FC(i), QFCA.FC(j)], ', ')};
QFCA.FC(j) = [];

end