function [ cp,d,m ] = closestpoint(F, T, mode)
% T target point set
% F floating point set
% For each point in F, find its closest point in T.
% i.e. F covers a smaller area.

if nargin<3
    mode = 1;
end

if size(F,2)~=size(F,2)
    error('Data size must be equal.')
end

if ~license('test', 'Statistics_Toolbox')
    mode = 1;
end

switch mode
    case 1
        kdOBJ = KDTreeSearcher(T);
        [m, d] = knnsearch(kdOBJ,F);
        cp = T(m,:);
    case 2
        [cp, d, m] = match_distance(F,T);
end

function [cp,d,m] = match_distance(F,T)
% T target point set
% F floating point set
[n,dim] = size(F);

cp=zeros(n,dim);
d=zeros(n,1);
m=zeros(n,1);

xt = T(:,1);
yt = T(:,2);

if dim==2 % 2d
    for i=1:n
        dist = sqrt((F(i,1)-xt).^2+(F(i,2)-yt).^2);
        [d(i),m(i)] = min(dist);
        cp(i,:) = [xt(m(i)) yt(m(i))];
    end
else % 3d
    zt = T(:,3);
    for i=1:n
        dist = sqrt((F(i,1)-xt).^2+(F(i,2)-yt).^2+(F(i,3)-zt).^2);
        [d(i),m(i)] = min(dist);
        cp(i,:) = [xt(m(i)) yt(m(i)) zt(m(i))];
    end
end