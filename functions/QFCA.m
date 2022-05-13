function model = QFCA(model,rnd)
% edited by Seirana
%
% QFCA computes the table of fully coupling relations and the list of blocked
% reactions for a metabolic network specified by its stoichiometric matrix
% and irreversible reactions and also returns the reduced metabolic network.
%
% USAGE:
%     fModel = QFCA(model)
%
% INPUTS:
%     model: the metabolic network with fields:
%         * .S -    the associated sparse stoichiometric matrix
%         * .rev -  the 0-1 indicator vector of the reversible reactions
%         * .rxns - the cell array of reaction abbreviations
%         * .mets - the cell array of metabolite abbreviations
%         * .lb -   the doulbe array of reaction flux lower bound
%         * .ub -   the doulbe array of reaction flux upper bound
%
% OUTPUT:
%     fModel: the reduced metabolic network with fields:
%         * .S -         the stoichiometric matrix(edited)
%         * .rev -       the 0-1 indicator vector of the reversible reactions(edited)
%         * .rxns -      cell array of reaction abbreviations(edited)
%         * .mets -      cell array of metabolite abbreviations(edited)
%         * .lb -        the doulbe array of reaction flux lower bound(edited)
%         * .ub -        the doulbe array of reaction flux upper bound(edited)
%         * .FC -        llist of fully coupled reaction
%         * .rxnNumber - the number of reactions in the original model
%         * .metNumber - the number of metabolites in the original model
%
% .. Authors:
%        - Mojtaba Tefagh, Stephen P. Boyd, 2019, Stanford University

r = rTOrnd(rnd);

model.rev = zeros(size(model.rxns,1),1);
for i = 1:size(model.rxns,1)
    if model.lb(i) < 0
        model.rev(i,1) = 1;
    end
end

QFCA.rev = model.rev;
QFCA.S = sparse(model.S);
QFCA.lb = model.lb;
QFCA.ub = model.ub;
QFCA.rxnN = string(model.rxnNumber);
QFCA.rxns = string(model.rxns);
QFCA.FC = 1:size(QFCA.S,2);
QFCA.FC = string(QFCA.FC');

QFCA.fctable = zeros(size(QFCA.S,2));

%%% identifying the dead-end metabolites and removing them from the network
[QFCA.S, metNum, ~] = unique(QFCA.S, 'rows', 'stable');
QFCA.mets = model.mets(metNum);

% aggregating all the isozymes
[~, reacNum, duplicates] = unique([QFCA.S.', QFCA.rev], 'rows', 'stable');
QFCA.S = QFCA.S(:, reacNum);

rN = zeros(size(reacNum));
rN(reacNum) = 1;

for i = size(rN,1):-1:1
    if rN(i) == 0
        QFCA.rxns(duplicates(i)) = {strjoin([QFCA.rxns(duplicates(i)), QFCA.rxns(i)], ', ')};
        QFCA.rxns(i) = [];
        
        QFCA.rxnN(duplicates(i)) = {strjoin([QFCA.rxnN(duplicates(i)), QFCA.rxnN(i)], ', ')};
        QFCA.rxnN(i) = [];
        
        QFCA.FC(duplicates(i)) = {strjoin([QFCA.FC(duplicates(i)), QFCA.FC(i)], ', ')};
        QFCA.FC(i) = [];
        
    end
end

QFCA.rev = QFCA.rev(reacNum);
QFCA.lb = QFCA.lb(reacNum);
QFCA.ub = QFCA.ub(reacNum);

fullCouplings = reacNum(duplicates);

%%% identifying the fully coupled pairs of reactions
% finding the trivial full coupling relations
[m, n] = size(QFCA.S);
flag = true;
while flag
    flag = false;
    for i = m:-1:1
        nzcols = find(QFCA.S(i, :));
        
        % check to see if the i-th row of S has only 2 nonzero elements
        if length(nzcols) == 2
            n = n-1;
            if (QFCA.ub(nzcols(2)) == QFCA.lb(nzcols(2))  && QFCA.ub(nzcols(1)) == QFCA.lb(nzcols(1))) || ...
                    (QFCA.ub(nzcols(2)) ~= QFCA.lb(nzcols(2))  && QFCA.ub(nzcols(1)) ~= QFCA.lb(nzcols(1)))
                QFCA = mergeFullyCoupled(QFCA, nzcols(1), nzcols(2), -QFCA.S(i, nzcols(1))/QFCA.S(i, nzcols(2)));
                fullCouplings(fullCouplings == reacNum(nzcols(2))) = reacNum(nzcols(1));
                reacNum(nzcols(2)) = [];
                flag = true;
            end
        end
    end
end

% finding the rest of full coupling relations
[Q, R, P] = qr(QFCA.S.');
tol = norm(QFCA.S, 'fro')*eps(class(QFCA.S));

% Z is the kernel of the stoichiometric matrix
rankS  = sum(abs(diag(R)) > tol);
Z = Q(:, rankS+1:n);
X = tril(Z*Z.');
Y = diag(diag(X).^(-1/2));
X = Y*X*Y;
for i = n:-1:2
    
    % j is the candidate reaction to be fully coupled to reaction i
    [M, j] = max(abs(X(i, 1:i-1)));
    
    % this is in fact cauchy-schwarz inequality
    if M > 1 - tol
        if (QFCA.ub(i) == QFCA.lb(i)  && QFCA.ub(j) == QFCA.lb(j)) || ...
                (QFCA.ub(i) ~= QFCA.lb(i)  && QFCA.ub(j) ~= QFCA.lb(j))
            QFCA = mergeFullyCoupled(QFCA, j, i, sign(X(i, j))*Y(j, j)/Y(i, i));
            fullCouplings(fullCouplings == reacNum(i)) = reacNum(j);
            reacNum(i) = [];
        end
    end
end
QFCA.S(:, QFCA.rev == -1) = -QFCA.S(:, QFCA.rev == -1);
[p, ~, ~] = find(P);
QFCA.S = QFCA.S(p(1:rankS), :);
QFCA.mets = QFCA.mets(p(1:rankS));
QFCA.metN = p(1:rankS);

QFCA.FCref = zeros(size(QFCA.fctable,1),1);
for i = 1:size(QFCA.FC,1)
    FC = str2double(split(QFCA.FC(i),','));
    QFCA.fctable(FC,FC) = 1;
    QFCA.FCref(FC,1) = i;
end

QFCA.S = full(QFCA.S);
model.QFCA = QFCA;

%%% QFCA finds the flux coupling coefficients
[~, n] = size(model.S);
k = n;
reacs = 1:n;
reactions = false(n, 1);
A = zeros(n);

fc = zeros(n,1);
for i = 1:n
    if sum(QFCA.fctable(i,:)) > 1
        fc(i) = 1;
    end
end

for i = k:-1:1
    if fc(i) == 0
        result = directionallyCoupled(model.S, model.rev, i, model.lb, model.ub);%(model.S, rev, i, model.lb, model.ub);
        dcouplings = result ~= 0;
        dcouplings(fc == 1) = false;
        dcouplings(i) = false;
        if any(dcouplings)
            %%% Irev ---> Irev flux coupling relations
            A(reacs, i) = 3*dcouplings;
            A(i, reacs) = 4*dcouplings.';
            dcouplings(i) = true;
            % correcting the reversibility conditions
            model.rev(i) = 0;
            % inferring by the transitivity of directional coupling relations
            A(reacs(model.rev == 1), i) = max(A(reacs(model.rev == 1), ...
                reacs(dcouplings)), [], 2);
            A(i, reacs(model.rev == 1)) = max(A(reacs(dcouplings), ...
                reacs(model.rev == 1)), [], 1);
            
            %%% Prev ---> Irev flux coupling relations
            if any(A(reacs(model.rev == 1), i) == 0)
                numLE = numLE + 1;
                coupled = false(n, 1);
                [Q, R, ~] = qr(transpose(model.S(:, ~dcouplings)));
                tol = norm(model.S(:, ~dcouplings), 'fro')*eps(class(model.S));
                Z = Q(model.rev(~dcouplings) == 1 & A(reacs(~dcouplings), i) == 0, ...
                    sum(abs(diag(R)) > tol)+1:end);
                coupled(~dcouplings & model.rev == 1 & A(reacs, i) == 0) = diag(Z*Z.') < tol^2;
                A(reacs(coupled), i) = 3;
                A(i, reacs(coupled)) = 4;
                % -1 indicates an uncoupled pair for remembering to skip
                % it without any need for further double check later
                A(reacs(~coupled & model.rev == 1 & A(reacs, i) == 0), ...
                    reacs(dcouplings)) = -1;
            end
            
            %%% metabolic network reductions induced by DCE
            if result(i) < 0
                model.S(:, i) = -model.S(:, i);
            end
            reactions(i) = true;
            for j = i+1:k
                if reactions(j)
                    if all(A(i, reacs(model.rev == 0 & ~reactions(reacs))) == ...
                            A(j, reacs(model.rev == 0 & ~reactions(reacs))))
                        A(i, j) = 2;
                        A(j, i) = 2;
                    elseif all(A(i, reacs(model.rev == 0 & ~reactions(reacs))) ...
                            <= A(j, reacs(model.rev == 0 & ~reactions(reacs))))
                        A(i, j) = 3;
                        A(j, i) = 4;
                    elseif all(A(i, reacs(model.rev == 0 & ~reactions(reacs))) ...
                            >= A(j, reacs(model.rev == 0 & ~reactions(reacs))))
                        A(i, j) = 4;
                        A(j, i) = 3;
                    end
                end
            end
        end
    end
end

% the usage of -1 was temporary and we return to our earlier convention
A(A == -1) = 0;
A(logical(eye(k))) = 1;

%%% postprocessing to fill in the flux coupling table for the original
% metabolic network from the flux coupling relations for the reduced one
a = 1:k;
map = repmat(fullCouplings.', k, 1) == repmat(a', 1, length(duplicates)); %reacNum
fctable = map.'*A*map;

fctable(logical(eye(size(fctable)))) = 1;

model.QFCA.DC_fctable = fctable;
end