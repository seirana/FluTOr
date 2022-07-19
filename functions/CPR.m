function model = CPR(mdl1, rnd, file_name)
% CPR find the coupling relation between each pair of reactions with F2C2 and QFCA methods
%
% USAGE:
%     model = CPR(mdl1, rnd, file_name)
%
% INPUT:
%     model: the metabolic network with fields:
%         * .S -    the associated sparse stoichiometric matrix
%         * .rev -  the 0-1 indicator vector of the reversible reactions
%         * .rxns - the cell array of reaction abbreviations
%         * .mets - the cell array of metabolite abbreviations
%         * .lb -   the doulbe array of reaction flux lower bound
%         * .ub -   the doulbe array of reaction flux upper bound
%
% OUTPUT:
%     model: the metabolic network with fields:
%         * .S -    the associated sparse stoichiometric matrix
%         * .rev -  the 0-1 indicator vector of the reversible reactions
%         * .rxns - the cell array of reaction abbreviations
%         * .mets - the cell array of metabolite abbreviations
%         * .lb -   the doulbe array of reaction flux lower bound
%         * .ub -   the doulbe array of reaction flux upper bound
%     model.QFCA: the reduced metabolic network with fields:
%         * .S -         the stoichiometric matrix(edited)
%         * .rev -       the 0-1 indicator vector of the reversible reactions(edited)
%         * .rxns -      cell array of reaction abbreviations(edited)
%         * .mets -      cell array of metabolite abbreviations(edited)
%         * .lb -        the doulbe array of reaction flux lower bound(edited)
%         * .ub -        the doulbe array of reaction flux upper bound(edited)
%         * .FC -        llist of fully coupled reaction
%         * .rxnNumber - the number of reactions in the original model
%         * .metNumber - the number of metabolites in the original model   


model = QFCA(model, rnd); % modified QFCA functions
model.F2C2.fctable = FC_F2C2(model);
model = QFCA_F2C2_subscription(model);

f_name = strcat(file_name, '_lb=', string(lb_bio));
save(strcat(adrs, f_name, '_model.mat'),'model'); % save model
end
