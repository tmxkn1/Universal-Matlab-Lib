function [ cp,d,m ] = closestpoint(F, T, mode)
%CLOSESTPOINT calculates the closest points between two given point sets
%
% C = closestpoint(F, T). 
%
% C = closestpoint(F, T, M) calculates the closest point C in T for each 
% point in F using one of the two methods:
%       1. kdtree method from statistics toolbox (M = 1)
%       2. brutal force (M = 2)
%
% To improve efficiency, make sure T is the bigger point set.
%
% C = closestpoint(F, T) omitting M to allow the function choose the 
% supposedly most efficient method automatically.
%
% [C, D] = closestpoint(...) also returns the distance between each point
% and corresponding closest point.
%
% [C, D, M] = closestpoint(...) also returns the indices of C in T.
%
% Copyright 2018, Zhengyi Jiang

if nargin < 2
    error('Not enough input arguments.');
end

if size(T,2)~=size(F,2)
    error('closestpoint(F, T) F and T must have the same number of columns.')
end

if nargin<3
    mode = 1;
    
    if size(T,1) < 20
        % point size is too small to make a good use of kd tree.
        mode = 2;
    end
end

if ~license('test', 'Statistics_Toolbox')
    mode = 2;
end

switch mode
    case 1
        kdOBJ = KDTreeSearcher(T);
        [m, d] = knnsearch(kdOBJ,F);
        cp = T(m,:);
    case 2
        [cp, d, m] = bfsearch(F,T);
end

%% brutal force
function [cp,d,m] = bfsearch(F,T)
% T target point set
% F floating point set
[n,dim] = size(F);

cp=zeros(n,dim);
d=zeros(n,1);
m=zeros(n,1);

for i = 1:n
    dist = 0;
    for j = 1:dim
        dist = dist + (F(i,j)-T(:,j)).^2;
    end
    [d(i),m(i)] = min(sqrt(dist));
    cp(i,:) = T(m(i),:);
end