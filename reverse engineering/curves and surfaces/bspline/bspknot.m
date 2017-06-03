function [knotV] = bspknot(ctrlPNum, k, segId, param)
if nargin < 3
    isUniform = true;
else
    isUniform = false;
end
n = ctrlPNum - 1;

if isUniform
    % Define uniform knot vector. (clapmed)
    knotV = zeros(1,n+k+1);
    knotV(k+1:n+1) = 1:1:n-k+1;
    knotV(n+2:end) = n-k+2;
    knotV = knotV(:);
else
    knotV = zeros(1,n+k+1);
    knotV(k+1:n+2) = param(segId(2:end-1)-1);
    knotV(n+2:end) = max(param);
end
