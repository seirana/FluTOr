function fctable = FC_F2C2(model)
% FC_table builds a model with desired fields that F2C2 function needs then calls the F2C2 function
%
% USAGE:
%     fctable = F2C2_table(model)
%
% INPUTS:
%     model: the metabolic network with fields:
%         * .S -    the associated sparse stoichiometric matrix
%         * .rxns - the cell array of reaction abbreviations
%         * .mets - the cell array of metabolite abbreviations
%         * .lb -   the doulbe array of reaction flux lower bound
%         * .ub -   the doulbe array of reaction flux upper bound
%
% OUTPUT:
%     fctable: the table to show the realation between each two reactions; fully coupled, partailly coupled, directionally coupled, or uncoupled

% add F2C2 fields for the model
fc_mdl = model;
fc_mdl.stoichiometricMatrix = model.S;
fc_mdl.Reactions = model.rxns;
fc_mdl.Metabolites = model.mets;
fc_mdl.reversibilityVector = zeros(size(model.rxns,1),1);
fc_mdl.lb = model.lb;
fc_mdl.ub = model.ub;

for i = 1:size(fc_mdl.Reactions,1)
    if fc_mdl.lb(i) < 0 && fc_mdl.ub(i) > 0
        fc_mdl.reversibilityVector(i,1) = 1;
    end
end

% F2C2
solver = 'glpk';%'linprog';

[fctable, ~] = F2C2(solver, fc_mdl);%, 10^(-1*rnd));

for i = 1:size(fctable,1)
    fctable(i,i) = 1;
    for j = 1:size(fctable,2)
        if fctable(i,j) > 1
            fctable(i,j) = 0;
        end
    end
end
end
