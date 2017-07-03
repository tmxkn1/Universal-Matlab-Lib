function [ P_h ] = FuncICP_toHomo( P,o )
%cahange2homocord changes coordinates to homogeneouse coordinates
%   Detailed explanation goes here
[N,M]=size(P);
if o==1
    P_h=[P;ones(1,M)];
    
elseif o==2
    P_h=[P ones(N,1)];
    
else
    P_h=[P;ones(1,M)];
end

