function model = FBA_FVA(model, rnd, bio, lb_bio, Ain, Bin)

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
