function model = FBA_FVA(model, rnd, bio, lb_bio, Ain, Bin)
% FBA_FVA calculates the flux range s for the reactions, 
%     and remove blocked reactions and deadend metabolites fro mthe model
%
% USAGE:
%     model = FBA_FVA(model, rnd, bio, lb_bio, Ain, Bin)
%
% INPUTS:
%     model: the metabolic network with fields:
%         * .S -    the associated sparse stoichiometric matrix
%         * .c -    a binary vectore for objective function
%         * .rxns - the cell array of reaction abbreviations
%         * .mets - the cell array of metabolite abbreviations
%         * .lb -   the doulbe array of reaction flux lower bound
%         * .ub -   the doulbe array of reaction flux upper bound
%     rnd: flux accuracy equals (default: 1e-5)
%     bio: biomass reaction number
%     lb_bio: lowest possible flux value for biomass reaction in the model
%     Ain: inequality matrix of linprog function(m*n)
%     Bin: right side of inequality matrix of linprog function(m*1)  Ain <= Bin 
%          (check https://de.mathworks.com/help/optim/ug/linprog.html#d123e112212 for more information)
%
% OUTPUT:
%     model: the metabolic network with fields:
%         * .S -    the associated sparse stoichiometric matrix
%         * .c -    a binary vectore for objective function
%         * .rxns - the cell array of reaction abbreviations
%         * .mets - the cell array of metabolite abbreviations
%         * .lb -   the doulbe array of reaction flux lower bound
%         * .ub -   the doulbe array of reaction flux upper bound
%         * .rxnNumner: to show the reaction number on the original model after converintg the model to irreversible, 
%                       and before deleting the blocked reactions and deadend metabolites

r = rTOrnd(rnd);

model.c = zeros(size(model.c));
model.c(bio) = 1;

options = optimset('linprog');
options.Display = 'off';

%FBA
% calculating optimum biomass rate
A = [];
b = [];

if nargin > 4
    A = [A;Ain];
    b = [b;Bin];
end

Aeq = model.S;
beq = zeros(size(Aeq,1),1);

lb = model.lb;
ub = model.ub;

if lb_bio > 0
    obj = zeros(size(model.c));
    obj(model.c == 1)=-1; % biomass
    
    [~,Sol.f,Sol.stat,~] = linprog(obj,A,b,Aeq,beq,lb,ub,options);
    
    if Sol.stat~=1
        ub(model.c == 1,1)=0;
    else
        ub(model.c == 1,1)=Sol.f*-1;
    end
    
    lb(model.c == 1,1) = lb_bio*ub(model.c == 1,1);
end

% FVA
% calculating flux range for all reactions
for i = 1:size(model.S,2)
    obj = zeros(size(model.c));
    obj(i)=-1;
    
    [~,Sol.f,Sol.stat,~] = linprog(obj,A,b,Aeq,beq,lb,ub,options);
    if Sol.stat~=1
        model.ub(i,1)=0;
    else
        model.ub(i,1)=Sol.f*-1;
    end
    [~,Sol.f,Sol.stat,~] = linprog(-obj,A,b,Aeq,beq,lb,ub,options);
    if Sol.stat~=1
        model.lb(i,1)=0;
    else
        model.lb(i,1)=Sol.f;
    end
end

model.lb = round(model.lb,r);
model.ub = round(model.ub,r);

% make and iireversible model
model = convertToIrreversible(model);

% make *.rxnNumber field
model.rxnNumber = zeros(size(model.rxns,1),1);
for i = 1:size(model.rxns,1)
    model.rxnNumber(i,1) = i;
end

% make *.metNumber field
model.metNumber = zeros(size(model.mets,1),1);
for i = 1:size(model.mets,1)
    model.metNumber(i,1) = i;
end

% remove blocked reactions and metabolites
model = remove_blocked_reactions_metabolites(model);
end
