function [model, biomass] = model_ecoliijo1366_davidi(file_name, rnd, condition)
% model_ecoliijo1366_davidi reads the original model and modifies based on the conditions of the problem
%
% USAGE:
%     [model, biomass] = model_ecoliijo1366_davidi(file_name, rnd, condition)
%
% INPUTS:
%     file_name: the name of the original model to load
%     rnd: flux accuracy (default: 1e-5)
%     condition: the number of the sub-model


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
carbon_source = xlsread(strcat(file_name, '_CarbonSources.xlsx'));

% block all carbon sources / desired one
f = find(string(model.rxns) == condition(1,4));
for t = 1:size(carbon_source,1)
    if carbon_source(t,1) ~= f
        model.lb(carbon_source(t,1)) = 0;
    else
        % active one carbon source
        if ~ismissing(condition(1,5))
            model.lb(f) = -1*str2double(condition(1,5));
            model.ub(f) = -1*str2double(condition(1,5));
        else
            model.lb(f) = -10;
        end
    end
end

% active biomass
bio = find(string(model.rxns) == 'BIOMASS_Ec_iJO1366_WT_53p95M');
biomass = model.rxns(bio);

% block other biomass rxns
model.lb(string(model.rxns) == 'BIOMASS_Ec_iJO1366_core_53p95M') = 0;
model.ub(string(model.rxns) == 'BIOMASS_Ec_iJO1366_core_53p95M') = 0;

% growth rate
model.lb(bio) = str2double(condition(1,3));
model.ub(bio) = str2double(condition(1,3));

model = FBA_FVA(model, rnd, bio, []);
end
