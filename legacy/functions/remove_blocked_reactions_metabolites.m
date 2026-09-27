function model = remove_blocked_reactions_metabolites(model)
% remove_blocked_reactions_metabolites removes blocked reactions and deadend metabolites
%
% USAGE:
%     model = remove_blocked_reactions_metabolites(model)
%
% INPUTS:
%     model: the metabolic network with fields:
%         * .lb - the lower bound of reaction fluxes
%         * .ub - the upper bound of reaction fluxes
%         * .rxns - reaction names
%         * .mets - metabolite names
%         * .S - stoichiometric matrix
%         * .c - objective coefficients
%         * .rxnNames - compelete reaction names
%         * .metNames - compelete reaction names
%         * .subSystems - subsystem which reaction belongs to
%         * .rxnNumber- reaction number
%         * .metNumber - metabolite number
%
% OUTPUT:
%     model: the metabolic network with fields:
%         * .lb - the lower bound of reaction fluxes
%         * .ub - the upper bound of reaction fluxes
%         * .rxns - reaction names
%         * .mets - metabolite names
%         * .S - stoichiometric matrix
%         * .c - objective coefficients
%         * .rxnNames - compelete reaction names
%         * .metNames - compelete reaction names
%         * .subSystems - subsystem which reaction belongs to
%         * .rxnNumber- reaction number
%         * .metNumber - metabolite number
%
% change negative irreversible rxns to positive
for i = 1:size(model.rxns,1)
    if model.lb(i) < 0 && model.ub(i) <= 0
        tmp = model.lb(i);
        model.lb(i) = -1*model.ub(i);
        model.ub(i) = -1*tmp;
        model.S(:,i) = -1*model.S(:,i);
    end
end

%%% remove blocked reactions from model
sz = size(model.rxns,1);
for i = sz:-1:1
    sm = 0;
    for j = 1:size(model.mets,1)
        sm = sm + abs(model.S(j,i));
    end
    if sm == 0 || (model.lb(i) == 0 && model.ub(i) == 0)
        model.S(:,i) = [];
        model.rxns(i) = [];
        model.lb(i) = [];
        model.ub(i) = [];
        model.c(i) = [];
        model.rxnNames(i) = [];
        model.subSystems(i) = [];
        model.rxnNumber(i) = [];
    end
end

%%% remove blocked metabolites from model
if size(model.rxns,2) > 0
    sz = size(model.mets,1);
    for i = sz:-1:1
        sm = 0;
        for j = 1:size(model.rxns,1)
            sm = sm + abs(model.S(i,j));
        end
        if sm == 0
            model.S(i,:) = [];
            model.mets(i) = [];
            model.b(i) = [];
            model.metNames(i) = [];
            model.metNumber(i) = [];
        end
    end
end
end
