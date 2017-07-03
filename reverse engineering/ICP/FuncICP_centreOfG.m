function [mu]=FuncICP_centreOfG(P)
%   This function calculates the center of mass for all points in P
%   P should be N x M matrix with M points in N-th dimension.

[N,M]=size(P);

mu=zeros(N,1);

mu=sum(P,2);
mu=mu/M;