function r = rTOrnd(rnd) 
% rTOrnd read flux accuracy value(rnd) and return the number of digits after decimal point
%
% USAGE:
%     r = rTOrnd(rnd) 
%
% INPUTS:
%     rnd: the flux accuracy
%
% OUTPUT:
%     r: the number of digits after zero

r = 0;
while rnd < 1
    rnd = rnd*10;
    r = r+1;
end
r = r+1;
end
