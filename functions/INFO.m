function [rxns, irr_lst, fix_lst, bio_lst, bio_QFCA, bio] = INFO(model, biomass)
% INFO reads the model and prepars the list of fixed and irreversible reactions
%
% USAGE:
%     [rxns, irr_lst, fix_lst, bio_lst, bio_QFCA, bio] = INFO(model, biomass)
%
% INPUTS:
%     model: the metabolic network with fields:
%         * .QFCA.rxns - the cell array of list of fully coupled reactions
%         * .rxns - the cell array of reaction abbreviations
%         * .QFCA.FC - a table that shows each two reactions are fully coupled or not
%         * .lb -   the doulbe array of reaction flux lower bound
%         * .ub -   the doulbe array of reaction flux upper bound
%     biomass: the name of biomass reaction
%
% OUTPUT:
%     rxns: a cell array that the first column refers to reaction numbers adn the second column refers to reaction names 
%     irr_lst: a list of the irreversible reactions with variable flux
%     fix_lst: a list of the reactions with  a fixed fluxv value
%     bio_lst: the list of reactions that fully coupled to biomass reaction
%     bio_QFCA: the index that shows this column of model.QFCA.rxns belng to biomass reaction and its fully couples reactions
%     bio: the number of the reaction 

bio_QFCA = strfind(string(model.QFCA.rxns),biomass);
if iscell(bio_QFCA)
    for i = 1:size(bio_QFCA,1)
        if ~isempty(bio_QFCA{i,1})
            bio_QFCA = i;
            break;
        end
    end
end

bio_lst = zeros(0,0);
l = model.QFCA.FC(bio_QFCA);
l = str2double(split(l,','));
bio_lst = [bio_lst;l];
bio_lst = sort(bio_lst);

% list of reactions with fixed flux
fix_lst = setdiff(find(model.lb == model.ub),bio_lst);
fix_lst = sort(fix_lst);

% list of reactions with irreversilbe varialbe flux
irr_lst = setdiff(find(model.lb ~= model.ub),bio_lst);
irr_lst = sort(irr_lst);

rxns = strings(size(model.rxns,1)+1,2);
rxns(end,1) = "trade_off";
rxns(end,2) = "length";

rxns(1:end-1,1) = model.rxnNumber;
rxns(1:end-1,2) = model.rxns;

bio = find(string(model.rxns) == biomass);
end
