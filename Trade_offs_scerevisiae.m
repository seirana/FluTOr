clc;clear;
% read A. thaliana model and find trade-offs
file_name = 'yeastGEM';
rnd = 1e-5; % flux accuracy equals 1e-5

[carbon_source,cb] = xlsread(strcat(adr, file_name, '_CarbonSources.xlsx')); % list of carbon sources
cb = string(cb(2:end,3));

% find trade-offs for each active carbon source
for c = 1:size(carbon_source,1)  % loop for active carbon source
    for j = 1:3 % loop for biomass lower bound
        lb_bio = 0.85+j*0.05; % biomass lower bound        
        [model, biomass] = model_scerevisiae_carbonsource(adr, file_name, rnd, lb_bio, carbon_source, c); % read the model
        model = CPR(model, rnd, file_name)  % find coupling relation between each pair of reactions          
        tradeoff_seaker(adrs, f_name, model, biomass); % find_tradeoffs
    end
end
