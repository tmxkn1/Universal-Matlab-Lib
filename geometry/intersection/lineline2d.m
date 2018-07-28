function [I, C]=lineline2d(l1p, l1n, l2p, l2n)
%oUTPUTS:
%
% I is the point of intersection 
%
% C indicates intersection condition: 
%   1 => parallel or colinear
%   2 => intersects

t = (l2p - l1p) / [l1n; -l2n];

if any(abs(t) == [inf inf])
    C = 1;
    I = nan(1, 2);
    return;
end

C = 2;
I = t(1) * l1n + l1p;