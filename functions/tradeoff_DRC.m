function [match, lst] = tradeoff_DRC(QFCA, F2C2, trade_offs, names)

s = size(trade_offs,1)-1;
DC = zeros(s,3);
trd = str2double(trade_offs(1:end-1,3:end));
rxns = trade_offs(1:end-1,2);

for i = 1:s
    if min(trd(i,:)) < 0
        DC(i,1) = 1;
    end
end
DC(:,2) = QFCA;
DC(:,3) = F2C2;

stat = zeros(s,1);
for i = 1:s
    stat(i,1) = 25*DC(i,1)+5*DC(i,2)+DC(i,3);
end
[GC,GR] = groupcounts(stat);
GC = 100*GC/s;
lst = zeros(50,1);
lst(GR+1,1) = GC;

for k = 1:size(rxns,1)
    if endsWith(rxns(k,1),'_f') || endsWith(rxns(k,1),'_b') || endsWith(rxns(k,1),'_r')
        aq = char(rxns(k,1));
        aq = aq(1,1:end-2);
        rxns(k,1) = aq;
    end
end

match = zeros(size(names,1),2);
for i = 1:s
    f = find(names == rxns(i,1));
    if match(f,1) < 1
        if min(trd(i,:)) < 0
            if (DC(i,2) == 0 || DC(i,2) == 4) &&  match(f,1) == 0
                match(f,1) = -1;
            end
            if DC(i,2) == 2 || DC(i,2) == 3
                match(f,1) = 1;
            end
            if (DC(i,3) == 0 || DC(i,3) == 4) &&  match(f,1) == 0
                match(f,2) = -1;
            end
            if DC(i,3) == 2 || DC(i,3) == 3
                match(f,2) = 1;
            end
        end
        if min(trd(i,:)) == 0
            if (DC(i,2) == 2 || DC(i,2) == 3) &&  match(f,1) == 0
                match(f,1) = -1;
            end
            if (DC(i,3) == 2 || DC(i,3) == 3) &&  match(f,1) == 0
                match(f,2) = -1;
            end
        end
    end
end

end