function [model, biomass] = model_ecoliijo1366_carbonsource(file_name, rnd, lb_bio, carbon_source, c)
%model_ecoliijo1366_carbonsource reads the original model and modifies based on the conditions of the problem
%
% USAGE:
%     [model, biomass] = model_ecoliijo1366_carbonsource(file_name, rnd, lb_bio, carbon_source, c)
%
% INPUTS:
%     file_name: the name of the original model to load
%     rnd: flux accuracy (default: 1e-5)
%     lb_bio: the rio of lower bound of lower bound of biomass to the optimum biomass flux value
%     carbon_source: the list of carbon sources
%     c: the active carbon source

% OUTPUT:
%     model: the metabolic network with fields:
%         * .lb - the lower bound of reaction fluxes
%         * .ub - the upper bound of reaction fluxes
%         * .rxns - reaction names
%         * .mets - metabolite names
%         * .S - stoichiometric matrix
%         * .c - objective coefficients
%     biomass: the name of active biomass reaction

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
