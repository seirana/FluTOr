function model = roundSmatrix(model, which)
if which == "org"
    S = model.S;
    S(S == 0) = Inf;    
    mini = min(abs(S));
    mm = min(mini);
    r = 1;
    while mm < 1
        mm = mm*10;
        r = r*10;
    end
    model.S(:,mini < 1) = model.S(:,mini < 1)*r;
    model.lb(mini < 0.001) = model.lb(mini < 0.001)/r;
    model.ub(mini < 0.001) = model.ub(mini < 0.001)/r;  

end
if which == "QFAC"
    model.QFAC.S(abs(model.QFAC.S) > 1) = round(model.QFAC.S(abs(model.QFAC.S) > 1),0);
    model.QFAC.S(0 < model.QFAC.S & model.QFAC.S < 0.001) = 0.001;
    model.QFAC.S(-0.001 < model.QFAC.S & model.QFAC.S < 0) = -0.001;
    model.QFAC.S(abs(model.QFAC.S) < 1) = round(model.QFAC.S(abs(model.QFAC.S) < 1),3);
end
end
