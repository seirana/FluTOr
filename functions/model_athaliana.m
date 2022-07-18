function [model, biomass] = model_athaliana(file_name, rnd, lb_bio, bio, nit_amu)
% model_athaliana2 reads the original model and modifies based on the conditions of the problem
%
% USAGE:
%     [model, biomass] = model_athaliana(file_name, rnd, lb_bio, bio, nit_amu)
%
% INPUTS:
%     file_name: the name of the original model to load
%     rnd: flux accuracy (default: 1e-5)
%     lb_bio: the rio of lower bound of lower bound of biomass to the optimum biomass flux value
%     bio: the name of active biomass raction
%     nit_amu: the uptake ratio of nitrate to ammunium

% OUTPUT:
%     model: the metabolic network with fields:
%         * .lb - the lower bound of reaction fluxes
%         * .ub - the upper bound of reaction fluxes
%         * .rxns - reaction names
%         * .mets - metabolite names
%         * .S - stoichiometric matrix
%         * .c - objective coefficients

load(strcat(file_name,'.mat')); % load the model

Ain = zeros(2,size(model.S,2));
Bin = zeros(2,1);

% set nitrogen source
if string(nit_amu) == 'NO3'
    model.lb(string(model.rxns) == 'Im_NH4') = 0;
    model.ub(string(model.rxns) == 'Im_NH4') = 0;
    model.S(:, string(model.rxns) == 'Im_NH4') = 0;
end

if string(nit_amu) == 'NH4'
    model.lb(string(model.rxns) == 'Im_NO3') = 0;
    model.ub(string(model.rxns) == 'Im_NO3') = 0;
    model.S(:, string(model.rxns) == 'Im_NO3') = 0;
end

if string(nit_amu) == '50:50'
    
    Ain(1,string(model.rxns) == 'Im_NO3') = 1;
    Ain(1,string(model.rxns) == 'Im_NH4') = -1;
    Bin(1,1) = rnd;
    
    Ain(2,string(model.rxns) == 'Im_NO3') = -1;
    Ain(2,string(model.rxns) == 'Im_NH4') = 1;
    Bin(2,1) = rnd;
end

bio = find(string(model.rxns) == string(bio));
bios = 547:549;

% block other biomass rxns
model.lb(bios(bios ~= bio),1) = 0;
model.ub(bios(bios ~= bio),1) = 0;

biomass = model.rxns(bio);

model = FBA_FVA(model, rnd, bio, lb_bio, Ain, Bin);

end
