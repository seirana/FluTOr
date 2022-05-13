clc;clear;

file_name = 'yeastGEM';
rnd = 1e-5;

[carbon_source,cb] = xlsread(strcat(adr, file_name, '_CarbonSources.xlsx'));
cb = string(cb(2:end,3));
% carbonsource
for c = 1:size(carbon_source,1)
    for j = 1:3
        lb_bio = 0.85+j*0.05;
        
        [mdl1, biomass] = model_scerevisiae_carbonsource(adr, file_name, rnd, lb_bio, carbon_source, c); % read the model
        
        mdl2 = QFCA(mdl1, rnd); % modified QFCA functions
        mdl2.F2C2.fctable = FC_F2C2(mdl2);
        model = QFCA_F2C2_subscription(mdl2);
        
        f_name = strcat(file_name, '_CB=', string(cb(c,1)), '_bioLB=', string(lb_bio));
        save(strcat(adrs, f_name, '_model.mat'),'model'); % save model
        
        tradeoff_seaker(adrs, f_name, model, biomass); % find_tradeoffs
    end
end
