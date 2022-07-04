function [model, biomass] = model_ecoliijo1366_carbonsource(file_name, rnd, lb_bio, carbon_source, c)

load(strcat(file_name,'.mat')); % load the model

% block all carbon sources
for t = 1:size(carbon_source,1)
    model.lb(carbon_source(t,1)) = 0;
end
% active one carbon source
model.lb(carbon_source(c,1)) = -10;

% active biomass
bio = find(string(model.rxns) == 'BIOMASS_Ec_iJO1366_WT_53p95M');
biomass = model.rxns(bio);

% block other biomass rxns
model.lb(string(model.rxns) == 'BIOMASS_Ec_iJO1366_core_53p95M') = 0;
model.ub(string(model.rxns) == 'BIOMASS_Ec_iJO1366_core_53p95M') = 0;

model = FBA_FVA(model, rnd, bio, lb_bio);
end
