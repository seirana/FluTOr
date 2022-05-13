function tradeoff_seaker(file_name, model, biomass)
if size(model.rxns,1) > 1
    
    [rxns, model, irr_lst, fix_lst, bio_lst, bio_QFCA, bio] = INFO(model, biomass);
    
    % Aeq = [k, a, u, s]
    % k*N = a; a_fix(no limitation); a_bio > 0; a_irr <= 0;
    trds = [];
    lim = 1000;
    M = lim+1;
    
    N = model.S;
    r = size(N,2);
    m = size(N,1);
    
    v = size(irr_lst,1);
    
    Aeq = [N', -1*eye(r), zeros(r,v), zeros(r,v)]; % k*N = a
    beq = zeros(r,1);
    
    Ain0 = [zeros(1,m), zeros(1,r), zeros(1,v), -ones(1,v)]; % 2 <= sum s = trade-off length
    bin0 = -2;
    
    cnd = zeros(1,r);
    cnd(1,bio_lst) = -1;
    Ain1 = [zeros(1,m), cnd, zeros(1,v), zeros(1,v)]; % 1 <= a_bio
    bin1 = -1;
    
    cnd = zeros(v,r);
    for i = 1:v
        cnd(i,irr_lst(i)) = 1;
    end
    
    Ain2 = [zeros(v,m), cnd, -eye(v), zeros(v)]; % a_irr <= u
    bin2 = zeros(v,1);
    
    Ain3 = [zeros(v,m), -cnd, -eye(v), zeros(v)]; % -u <= a_irr
    bin3 = zeros(v,1);
    
    Ain4 = [zeros(v,m), cnd, zeros(v), M*eye(v)]; % a_irr + M*s <= M-1
    bin4 = (M-1)*ones(v,1);
    
    Ain5 = [zeros(v,m), -cnd, zeros(v), -M*eye(v)]; % 0 <= a_irr_i + M*s
    bin5 = zeros(v,1);
    
    Ain = [Ain0;Ain1;Ain2;Ain3;Ain4;Ain5];
    bin = [bin0;bin1;bin2;bin3;bin4;bin5];
    
    for i = 1:size(model.QFCA.FC,1)
        if i ~= bio_QFCA % ~bio
            l = model.QFCA.FC(i);
            l = str2double(split(l,','));
            if model.lb(l(1)) ~= model.ub(l(1)) % ~fix
                nAin = zeros(1,size(Ain,2));
                nbin = zeros(1,size(bin,2));
                for j = 1:size(l,1)
                    f = find(irr_lst == l(j));
                    nAin(1,m+r+v+f) = 1;
                end
                nbin = 1;
                Ain = [Ain;nAin];
                bin = [bin;nbin];
            end
        end
    end
    
    a_lb = zeros(r,1);
    a_ub = zeros(r,1);
    
    a_lb(bio_lst,1) = 0;
    a_ub(bio_lst,1) = Inf;
    
    a_lb(fix_lst,1) = -Inf;
    a_ub(fix_lst,1) = Inf;
    
    a_lb(irr_lst,1) = -lim;
    a_ub(irr_lst,1) = 0;
    
    %     k          a     u           s
    lb = [-Inf(m,1); a_lb; zeros(v,1); zeros(v,1)];
    ub = [ Inf(m,1); a_ub; lim*ones(v,1); ones(v,1)];
    
    intcon = size(lb,1)-v+1:size(lb,1); % integer irreversibles
    obj = [zeros(m,1); zeros(r,1); zeros(v,1); ones(v,1)]; % factors for objective function
    options = optimoptions('intlinprog', 'IntegerPreprocess','none','LPPreprocess','none' ,'MaxTime', 3600);
    
    while 1
        sol = intlinprog(obj,intcon,Ain,bin,Aeq,beq,lb,ub,options);
        sol = round(sol,6);
        if ~isempty(sol) && sum(abs(sol)) > 0
            
            s = sum(round(sol(m+r+v+1:end))); % trade-off length
            
            % new trade-off
            new_trd = irr_lst(round(sol(m+r+v+1:end)) == 1);
            
            % if there is any fully coupled reactions to reactions in the trade-off, find them
            elements = {};
            for i = 1:size(new_trd,1)
                fc = find(model.QFCA.fctable(new_trd(i),:) == 1);
                elements = [elements {fc}];
            end
            
            % find all combination of tradeoff among fully coupled rxns
            combinations = cell(1, numel(elements));
            [combinations{:}] = ndgrid(elements{:});
            combinations = cellfun(@(x) x(:), combinations,'uniformoutput',false);
            result = [combinations{:}]; % list of new trade-offs
            
            trd = zeros(size(model.rxns,1)+1,size(result,1));
            % add new conditions to the program
            for i = 1:size(result,1)
                s_in = zeros(1,v);
                for j = 1:size(result,2)
                    s_in(1,irr_lst == result(i,j)) = 1;
                    trd(result(i,j),i) = -1;
                    trd(bio,i) = 1;
                    trd(end,i) = s;
                end
                new_bin = [zeros(1,m), zeros(1,r), zeros(1,v), s_in];
                Ain = [Ain; new_bin];
                bin = [bin; s-1];
            end
            
            trds = [trds trd];
            
            if size(trds,1) > 0 && size(trds,2) > 0
                trade_offs = [rxns, num2cell(trds)];
                save(strcat(file_name, '_tradeoffs.mat'),'trade_offs');
            end
        else
            break;
        end
    end
end
end