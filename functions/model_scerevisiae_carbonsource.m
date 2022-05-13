function [model, biomass] = model_scerevisiae_carbonsource(file_name, rnd, lb_bio, carbon_source, c)

load(strcat(file_name,'.mat'));

% carbon source
% deactive all carbon sources / glucose
for t = 1:size(carbon_source,1)
    model.lb(carbon_source(t,1)) = 0;
end
% active glocuse
model.lb(carbon_source(c,1)) = -10; %'r_1714' = 'D-glucose exchange'

% active biomass
bio = find(string(model.rxnNames) == 'growth');
biomass = model.rxns(bio);

model = FBA_FVA(model, rnd, bio, lb_bio);
end