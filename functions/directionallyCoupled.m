function result = directionallyCoupled(S, rev, i, lb, ub)
% directionallyCoupled finds all the directionally coupled reactions to i
%
% USAGE:
%
%    [certificate, result] = directionallyCoupled(S, rev, i, solver)
%
% INPUTS:
%    S:         the associated sparse stoichiometric matrix
%    rev:       the 0-1 vector with 1's corresponding to the reversible reactions
%    i:         the index of reaction to which others are directionally coupled
%    solver:    the LP solver to be used; the currently available options are
%               'gurobi', 'linprog', and otherwise the default COBRA LP solver
%
% OUTPUTS:
%    certificate:    the fictitious metabolite for the positive certificate;
%                    S.'*certificate will be the corresponding directional
%                    coupling equation
%    result:         the result returned by the LP solver; all the -1 entries
%                    are directionally coupled to reaction i and the other
%                    entries except i are zero
%
% .. Authors:
%       - Mojtaba Tefagh, Stephen P. Boyd, 2019, Stanford University

[~, n] = size(S);
irevIndex = [1:i-1, i+1:n];
irevIndex = irevIndex(rev([1:i-1, i+1:n]) == 0);

A = [];
b = [];

Aeq = S;
beq = zeros(size(Aeq,1),1);

obj = zeros(n,1);
obj(irevIndex) = -1;

options = optimset('linprog');
options.Display = 'off';

if lb(i) == 0
    lb(i) = 1e-5;
end

[result,Sol.f,Sol.stat,~] = linprog(obj,A,b,Aeq,beq,lb,ub,options);



% [m, n] = size(S);
% irevIndex = [m+1:m+i-1, m+i+1:m+n];
% irevIndex = irevIndex(rev([1:i-1, i+1:n]) == 0);
% 
% model.A = [S.', -speye(n)];
% model.rhs = zeros(n, 1);
% 
% model.lb = [-Inf(m,1); lb];
% model.ub = [Inf(m,1); ub];
% 
% model.lb(m+i) = 0;
% model.ub(m+i) = 0;
% 
% problem.f = model.obj;
% problem.Aineq = model.A(rev == 0, :);
% problem.bineq = model.rhs(rev == 0);
% problem.Aeq = model.A(rev ~= 0, :);
% problem.beq = model.rhs(rev ~= 0);
% problem.lb = model.lb;
% problem.ub = model.ub;
% problem.solver = 'linprog';
% problem.options = optimset('Display', 'off');
% result = linprog(problem);
% result = result(m+1:end,1);

end