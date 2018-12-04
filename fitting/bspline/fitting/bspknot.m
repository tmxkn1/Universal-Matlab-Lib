function [knotV] = bspknot(ctrlPNum, k, segId, param)
% Calculates knots for a B-spline curve fitting.
%
% -------------------------------------------------------------------------
% USE:
%
% V = bspknot(C, K) calculates a uniform set of knots for a B-spline curve
%     of order K with control points, C.
%
%   Input:  C - an m*3 matrix, containing the X-, Y- and Z-coordinates of 
%              the control points of the B-spline curve.
%           K - a scalar, representing the order of the B-spline curve.
%
%   Output: V - a 1*k vector containing the knots of the B-spline curve.
%
% V = bspknot(C, K, S, U) calculates a non-uniform set of knots, V, for a 
%     B-spline curve of order K. C contains the controls points of the 
%     curve; S contains the indices of the points at which segments end;
%     and U contains the parametric coordinates on the B-spline curve.
%
% -------------------------------------------------------------------------
% See also: bspbasis, bspchordlparam, bspcontrolpt, bspgenpoints, 
% bspcurvefit

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
