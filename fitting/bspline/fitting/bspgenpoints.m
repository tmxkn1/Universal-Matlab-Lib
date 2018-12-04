function fitp = bspgenpoints(u, ctrlp, K, knotV)
% Generates points on a B-spline curve.
%
% -------------------------------------------------------------------------
% USE:
%
% F = bspgenpoints(U, C, K, V) 
%
% Input:
%   U - a 1*n vector containing the parametric coordinates on the B-spline 
%       curve.
%   C - an m*3 matrix, containing the X-, Y- and Z-coordinates of the 
%       control points of the B-spline curve.
%   K - a scalar, representing the order of the B-spline curve.
%   V - a 1*k vector containing the knots of the B-spline curve.
%
% Output:
%   F - an n*3 matrix, containing the X-, Y- and Z-coordinates of the
%       points generated on the B-spline curve.
%
% -------------------------------------------------------------------------
% See also: bspbasis, bspchordlparam, bspcontrolpt, bspknot, bspcurvefit

U = bspbasis(u, numel(ctrlp(:,1)), K, knotV);

FX=U*ctrlp(:,1);
FY=U*ctrlp(:,2);
FZ=U*ctrlp(:,3);

fitp=[FX,FY,FZ];