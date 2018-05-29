function U = bspbasis(param, ctrlPNum, K, knotV)
% Input arguments:
%   param:
%       a row vector of n number of parameters, i.e. u/w coordinates.
%
%   ctrlPNum:
%       a scalar, equals the number of control points.
%
%   K:
%       a scalar, equals the order of the curve.
%
%   knotV:
%       a row vector of no_of_ctrl number of knot vectors, can be obtained
%       by [knotV, d1, d2] = sub_Bspline_knotv(no_of_ctrl, deg_of_curve).
%
% Output arguments:
%   U:
%       a no_of_ctrl x 3 matrix; each column is the x-, y- and
%       z-coordinates of all control points.

uNum = numel(param);
uMax = max(param);
index = 1:ctrlPNum;

% preallocating variables
U = zeros(numel(param),ctrlPNum);
N = zeros(ctrlPNum+1,K);

% Calculate the denominators for basis functions (k>2). -may be useful when
% the size of point data is substantial, so the calculation is not repeated.
d1 = zeros(ctrlPNum,K);
d2 = d1;
for m=2:K
    d1(:,m) = knotV(index+m-1) - knotV(index);
    d2(:,m) = knotV(index+m) - knotV(index+1);
end

% optimising for loop
knotV = knotV(:); 
knotv1 = knotV(1:ctrlPNum+1); 
knotv2 = knotV(2:ctrlPNum+2);
knotVSize = size(knotv1);
knotc1 = knotV(ctrlPNum+1); knotc2 = knotV(ctrlPNum+2);

for ui = 1:uNum
    % k = 1
    u_ = param(ui)*ones(knotVSize);
    NA = u_>=knotv1 & u_<knotv2;
    if param(ui) == uMax && knotc1 == uMax && knotc2 == uMax
        NA(1:ctrlPNum+1) = 0;
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