function model = CPR(mdl1, rnd, file_name)

mdl2 = ComplexTradeoff(mdl1); % find complexes

mdl3 = QFCA(mdl2, rnd); % modified QFCA functions
mdl3.F2C2.fctable = FC_F2C2(mdl3);
model = QFCA_F2C2_subscription(mdl3);

f_name = strcat(file_name, '_lb=', string(lb_bio));
save(strcat(adrs, f_name, '_model.mat'),'model'); % save model
end