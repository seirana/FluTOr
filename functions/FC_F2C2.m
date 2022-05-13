function fctable = FC_F2C2(model)

%%% add F2C2 fields for the model
fc_mdl = model;
fc_mdl.stoichiometricMatrix = model.S;
fc_mdl.Reactions = model.rxns;
fc_mdl.Metabolites = model.mets;
fc_mdl.reversibilityVector = zeros(size(model.rxns,1),1);
fc_mdl.lb = model.lb;
fc_mdl.ub = model.ub;

for i = 1:size(fc_mdl.Reactions,1)
    if fc_mdl.lb(i) < 0 && fc_mdl.ub(i) > 0
        fc_mdl.reversibilityVector(i,1) = 1;
    end
end

%%% F2C2
solver = 'glpk';%'linprog';

[fctable, ~] = F2C2(solver, fc_mdl);%, 10^(-1*rnd));

for i = 1:size(fctable,1)
    fctable(i,i) = 1;
    for j = 1:size(fctable,2)
        if fctable(i,j) > 1
            fctable(i,j) = 0;
        end
    end
end
end
