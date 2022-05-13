function [rxns, model, irr_lst, fix_lst, bio_lst, bio_QFCA, bio] = INFO(model, biomass)

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

fix_lst = setdiff(find(model.lb == model.ub),bio_lst);
fix_lst = sort(fix_lst);

irr_lst = setdiff(find(model.lb ~= model.ub),bio_lst);
irr_lst = sort(irr_lst);

rxns = strings(size(model.rxns,1)+1,2);
rxns(end,1) = "trade_off";
rxns(end,2) = "length";

rxns(1:end-1,1) = model.rxnNumber;
rxns(1:end-1,2) = model.rxns;

bio = find(string(model.rxns) == biomass);
end
