function U = bspbasis(param, ctrlPCount, K, knotV)
% Calculates B-spline basis functions
%
% -------------------------------------------------------------------------
% USE:
%
% U = bspbasis(P, C, K, V)
%
%   Input:  P - a 1*n vector containing the parameter coordinates, i.e. 
%               u/w coordinates.
%           C - a scalar, representing the number of control points.
%           K - a scalar, representing the order of the curve.
%           V - a 1*m vector containing the knot vectors.
%
%   Output: U - an n*C matrix containing the coefficients of the basis 
%               functions.
%
% -------------------------------------------------------------------------
% See also: bspgenpoints, bspcontrolpt, bspknot, bspchordlparam,
% bspcurvefit

uNum = numel(param);
uMax = max(param);
index = 1:ctrlPCount;

% preallocating variables
U = zeros(uNum,ctrlPCount);
N = zeros(ctrlPCount+1,K);

% Calculate the denominators for basis functions (k>2). -may be useful when
% the size of point data is substantial, so the calculation is not repeated.
d1 = zeros(ctrlPCount,K);
d2 = d1;
for m=2:K
    d1(:,m) = knotV(index+m-1) - knotV(index);
    d2(:,m) = knotV(index+m) - knotV(index+1);
end

% optimising for loop
knotV = knotV(:); 
knotv1 = knotV(1:ctrlPCount+1); 
knotv2 = knotV(2:ctrlPCount+2);
knotVSize = size(knotv1);
knotc1 = knotV(ctrlPCount+1); knotc2 = knotV(ctrlPCount+2);

for ui = 1:uNum
    % k = 1
    u_ = param(ui)*ones(knotVSize);
    NA = u_>=knotv1 & u_<knotv2;
    if param(ui) == uMax && knotc1 == uMax && knotc2 == uMax
        NA(1:ctrlPCount+1) = 0;
        NA(end-1,1) = 1;
    end
    N(:,1) = NA;
    % k > 2
    for k = 2:K
        p1 = (param(ui) - knotV(index)) ./ d1(:,k) .* N(index,k-1);
        p1(isnan(p1)) = 0;
        p2 = (knotV(index+k) - param(ui)) ./ d2(:,k) .* N(index+1,k-1);
        p2(isnan(p2)) = 0;
        N(index,k) = p1 + p2;
    end
    U(ui,:)=N(1:end-1,k);
    U(ui,:)=N(1:end-1,k);
end