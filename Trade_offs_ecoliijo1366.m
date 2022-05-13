clc;clear;
warning off;

addpath(genpath('/home/seirana/cobratoolbox'));
addpath(genpath('/home/seirana/Desktop/Trade-offsReactions/Potential Trade-offs'));
addpath(genpath('/home/seirana/Desktop/Trade-offsReactions/F2C2 20v0.95b'));

ad = '/home/seirana/Desktop/Trade-offsReactions/Potential Trade-offs/';
adr = strcat(ad,'models/');
addpath(genpath(strcat(ad)));

file_name = 'Ecoli_iJO1366';
adrs = strcat(ad,file_name,'/');

rnd = 1e-5;

[~,~,conditions] = xlsread(strcat(adr,file_name, '_Conditions.xlsx'));
% davidi
for c = 2:32
    [mdl1, biomass] = model_ecoliijo1366_davidi(adr, file_name, rnd, string(conditions(c,:))); % read the model

    if size(mdl1.rxns,2) > 0
        f = find(string(mdl1.rxns) == biomass);
        if mdl1.ub(f) > 0

            mdl2 = ComplexTradeoff(mdl1); % find complexes

            mdl3 = QFCA(mdl2, rnd); % modified QFCA functions
            mdl3.F2C2.fctable = FC_F2C2(mdl3);
            model = QFCA_F2C2_subscription(mdl3);

            f_name = strcat(file_name, '_NEW_condition=', string(conditions(c,1)));
            save(strcat(adrs, f_name, '_model.mat'),'model'); % save model

            tradeoff_seaker(adrs, f_name, model, biomass); % find_tradeoffs
        end
    end
end

% [carbon_source,cb] = xlsread(strcat(adr, file_name, '_CarbonSources.xlsx'));
% cb = string(cb(2:end,2));
% 
% % carbonsource
% for c = 23:23%size(carbon_source,1)
%     for j = 2:3
%         if j < 4
%             lb_bio = 0.80+j*0.05;
%         end
%         if j == 4
%             lb_bio = 0.79+j*0.05;
%         end
%         
%         [mdl1, biomass] = model_ecoliijo1366_carbonsource(adr, file_name, rnd, lb_bio, carbon_source, c); % read the model
%         disp(size(mdl1.S));
%         f = find(string(mdl1.rxns) == biomass);
%         if mdl1.ub(f) > 057:58
%             mdl2 = ComplexTradeoff(mdl1); % find complexes
%             
%             model = QFCA(mdl2, rnd); % modified QFCA functions
%             
%             f_name = strcat(file_name, '_CB=', cb(c,1), '_bioLB=', string(lb_bio));
%             save(strcat(adrs, f_name, '_model.mat'),'model'); % save model
%             
%             tradeoff_seaker(adrs, f_name, model, biomass); % find_tradeoffs
%         end
%     end
% end
