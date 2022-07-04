clc;clear;
% reads E. coli model and finds trade-offs
file_name = 'Ecoli_iJO1366';
rnd = 1e-5; % flux accuracy equals 1e-5

[~,~,conditions] = xlsread(strcat(adr,file_name, '_Conditions.xlsx')); % reads a set of constraitns applies on the model
% find trade-offs under each set of constraitns
for c = 2:size(conditions,1) % loop for set of constraitns
    [mdl1, biomass] = model_ecoliijo1366_davidi(file_name, rnd, string(conditions(c,:))); % read the model
    
    mdl2 = QFCA(mdl1, rnd); % modified QFCA functions
    mdl2.F2C2.fctable = FC_F2C2(mdl2); % modified F2C2 function
    model = QFCA_F2C2_subscription(mdl2); % find mutual fully coupled reactions
    
    f_name = strcat(file_name, '_condition=', string(conditions(c,1)));
    save(strcat(adrs, f_name, '_model.mat'),'model'); % save model
    
    tradeoff_seaker(f_name, model, biomass); % find_tradeoffs
end

[carbon_source,cb] = xlsread(strcat(adr, file_name, '_CarbonSources.xlsx')); % list of carbon sources
cb = string(cb(2:end,2));

% find trade-offs for each active carbon source
for c = 1:size(carbon_source,1) % loop for active carbon source
    for j = 1:3 % loop for biomass lower bound
        lb_bio = 0.85+j*0.05; % biomass lower bound
        
        [mdl1, biomass] = model_ecoliijo1366_carbonsource(file_name, rnd, lb_bio, carbon_source, c); % read the model
        
        mdl2 = QFCA(mdl1, rnd); % modified QFCA functions
        mdl2.F2C2.fctable = FC_F2C2(mdl2); % modified F2C2 function
        model = QFCA_F2C2_subscription(mdl2); % find mutual fully coupled reactions
        
        f_name = strcat(file_name, '_CB=', cb(c,1), '_bioLB=', string(lb_bio));
        save(strcat(adrs, f_name, '_model.mat'),'model'); % save model
        
        tradeoff_seaker(f_name, model, biomass); % find_tradeoffs
    end
end
