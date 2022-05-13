% clc;clear;
warning off;

addpath(genpath('/home/seirana/cobratoolbox'));
addpath(genpath('/home/seirana/Desktop/Trade-offsReactions/Potential Trade-offs'));
addpath(genpath('/home/seirana/Desktop/Trade-offsReactions/F2C2 20v0.95b'));

ad = '/home/seirana/Desktop/Trade-offsReactions/Potential Trade-offs/';
adr = strcat(ad,'models/');
addpath(genpath(strcat(ad)));

file_name = 'ArabidopsisCoreModel';
adrs = strcat(ad,file_name,'/');

rnd =1e-5;
nit_amu = {'NO3','NH4','50:50'};
active_bio = {'Bio_CLim','Bio_NLim','Bio_opt'};

for b = 3:3
    for na = 1:3
        for j = 5:5
            if j < 6
                 lb_bio = 0.75+j*0.05;
            end
            if j == 6
                lb_bio = 0.69+j*0.05;
            end
            
            [mdl1, biomass] = model_athaliana(adr, file_name, rnd, lb_bio, active_bio(b), nit_amu(na)); % read the model
            
            mdl2 = ComplexTradeoff(mdl1); % find complexes
            
            mdl3 = QFCA(mdl2, rnd); % modified QFCA functions
            mdl3.F2C2.fctable = FC_F2C2(mdl3);
            model = QFCA_F2C2_subscription(mdl3);
            
            f_name = strcat(file_name, '_bio=', string(active_bio(b)), '_', string(nit_amu(na)), '_lb=', string(lb_bio));
            save(strcat(adrs, f_name, '_model.mat'),'model'); % save model
            
            tradeoff_seaker(adrs, f_name, model, biomass); % find_tradeoffs
        end
    end
end
