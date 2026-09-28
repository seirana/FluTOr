clc;clear;
% read A. thaliana model and find trade-offs
file_name = 'ArabidopsisCoreModel';
rnd =1e-5;  % flux accuracy equals 1e-5
nit_amu = {'NO3','NH4','50:50'}; % nitrogen source
active_bio = {'Bio_CLim','Bio_NLim','Bio_opt'}; % active biomass 

for b = 1:3 % loop for active biomass
    for na = 1:3 % loop for nictogen source
        for j = 1:3 % loop for value of biomass lower bound
            lb_bio = 0.85+j*0.05; % biomass lower bound            
            [model, biomass] = model_athaliana(file_name, rnd, lb_bio, active_bio(b), nit_amu(na)); % read the model            
            model = CPR(model, rnd, file_name)  % find coupling relation between each pair of reactions          
            tradeoff_seaker(f_name, model, biomass); % find_tradeoffs
        end
    end
end
