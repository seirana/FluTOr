function model = QFCA_F2C2_subscription(model)
%QFCA_F2C2 finds mutual fully coupled reactions betwen the two method QFCA and F2C2
%
% USAGE:
%     model = QFCA_F2C2_subscription(model)
%
% INPUTS:
%     model: the metabolic network with fields:
%         * .QFCA.fctable - a matrix that shows reaction between each two reactions defined by QFCA method
%         * .F2C2.fctable - a matrix that shows reaction between each two reactions defined by F2C2 method
%         * .QFCA.FCref - returns reference number for fullycpoupled reactions
%         * .QFCA.rxns - the list reaction names which are fully coupled
%         * .QFCA.rxnN - the list reaction numbers which are fully coupled
%         * .QFCA.FC- the list reaction names which are fully coupled
%
% OUTPUT:
%     model: the metabolic network with fields:
%         * .QFCA.fctable - a matrix that shows reaction between each two reactions defined by QFCA method
%         * .F2C2.fctable - a matrix that shows reaction between each two reactions defined by F2C2 method
%         * .QFCA.FCref - returns reference number for fullycpoupled reactions
%         * .QFCA.rxns - the list reaction names which are fully coupled
%         * .QFCA.rxnN - the list reaction numbers which are fully coupled
%         * .QFCA.FC- the list reaction names which are fully coupled

QFCA04 = model.QFCA.DC_fctable;
QFCA01 = model.QFCA.DC_fctable;
QFCA01(QFCA01 > 1) = 0;

model.QFCA = rmfield(model.QFCA,{'DC_fctable','lb','ub','S','rev','mets','metN'}); % remove some fields

F2C204 = model.F2C2.fctable;
F2C201 = model.F2C2.fctable;
F2C201(F2C201 > 1) = 0;

diff = abs(QFCA01 - F2C201);

model.QFCA.FCref = zeros(size(QFCA04,1),1);
model.QFCA.rxns = string(model.rxns);
model.QFCA.rxnN = string(model.rxnNumber);
model.QFCA.FC = 1:size(model.S,2);
model.QFCA.FC = string(model.QFCA.FC');

% update modeol.QFCA 
for i = 1:size(QFCA04,1)
    for j = i+1:size(QFCA04,2)
        if diff(i,j) == 1
            if QFCA04(i,j) == 1
                QFCA04(i,j) = 2; % update fctable
                QFCA04(j,i) = 2;
            end
            if F2C204(i,j) == 1
                F2C204(i,j) = 2;
                F2C204(j,i) = 2;
            end
        end
    end
    f = find(QFCA04(i,:) == 1);
    for j = 1:size(f,2)
        model.QFCA.FCref(f(1,j)) = f(1,1); % update FCref
    end
end

for i = size(model.QFCA.FCref,1):-1:1
    % update rxns, rxnN, FC
    t = model.QFCA.FCref(i);
    if t < i
        model.QFCA.rxns(t) = {strjoin([model.QFCA.rxns(t), model.QFCA.rxns(i)], ', ')};
        model.QFCA.rxns(i,:) = [];
        
        model.QFCA.rxnN(t) = {strjoin([model.QFCA.rxnN(t), model.QFCA.rxnN(i)], ', ')};
        model.QFCA.rxnN(i,:) = [];
        
        model.QFCA.FC(t) = {strjoin([model.QFCA.FC(t), model.QFCA.FC(i)], ', ')};
        model.QFCA.FC(i,:) = [];
    end
end
model.QFCA.fctable = QFCA04;
model.F2C2.fctable = F2C204;

end
