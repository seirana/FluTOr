function fctable = FC_F2C2(model)
% FC_F2C2 calculates the flux range s for the reactions, 
%     and remove blocked reactions and deadend metabolites fro mthe model
%
% USAGE:
%     fctable = FC_F2C2(model)
%
% INPUTS:
%     model: the metabolic network with fields:
%         * .S -    the associated sparse stoichiometric matrix
%         * .c -    a binary vectore for objective function
%         * .rxns - the cell array of reaction abbreviations
%         * .mets - the cell array of metabolite abbreviations
%         * .lb -   the doulbe array of reaction flux lower bound
%         * .ub -   the doulbe array of reaction flux upper bound
%     rnd: flux accuracy equals (default: 1e-5)
%     bio: biomass reaction number
%     lb_bio: lowest possible flux value for biomass reaction in the model
%     Ain: inequality matrix of linprog function(m*n)
%     Bin: right side of inequality matrix of linprog function(m*1)  Ain <= Bin 
%          (check https://de.mathworks.com/help/optim/ug/linprog.html#d123e112212 for more information)
%
% OUTPUT:
%     model: the metabolic network with fields:
%         * .S -    the associated sparse stoichiometric matrix
%         * .c -    a binary vectore for objective function
%         * .rxns - the cell array of reaction abbreviations
%         * .mets - the cell array of metabolite abbreviations
%         * .lb -   the doulbe array of reaction flux lower bound
%         * .ub -   the doulbe array of reaction flux upper bound
%         * .rxnNumner: to show the reaction number on the original model after converintg the model to irreversible, 
%                       and before deleting the bloc

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
