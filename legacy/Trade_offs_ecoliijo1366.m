clc;clear;
% reads E. coli model and finds trade-offs
file_name = 'Ecoli_iJO1366';
rnd = 1e-5; % flux accuracy equals 1e-5

[~,~,conditions] = xlsread(strcat(adr,file_name, '_Conditions.xlsx')); % reads a set of constraitns applies on the model
% find trade-offs under each set of constraitns
for c = 2:size(conditions,1) % loop for set of constraitns
    [model, biomass] = model_ecoliijo1366_davidi(file_name, rnd, string(conditions(c,:))); % read the model    
     model = CPR(model, rnd, file_name)  % find coupling relation between each pair of reactions    
    tradeoff_seaker(f_name, model, biomass); % find_tradeoffs
end

[carbon_source,cb] = xlsread(strcat(adr, file_name, '_CarbonSources.xlsx')); % list of carbon sources
cb = string(cb(2:end,2));

% find trade-offs for each active carbon source
for c = 1:size(carbon_source,1) % loop for active carbon source
    for j = 1:3 % loop for biomass lower bound
        lb_bio = 0.85+j*0.05; % biomass lower bound        
        [model, biomass] = model_ecoliijo1366_carbonsource(file_name, rnd, lb_bio, carbon_source, c); % read the model            
        model = CPR(model, rnd, file_name)  % find coupling relation between each pair of reactions    
        tradeoff_seaker(f_name, model, biomass); % find_tradeoffs
    end
end
