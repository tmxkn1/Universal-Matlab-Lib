function ctrlp = bspcontrolpt(data, u, ctrlPNum, K, knotV)
%Calculates a set of control points for a 3D B-spline approximation.
%
% -------------------------------------------------------------------------
% USE:
%
% C = bspcontrolpt(D, U, N, K, V)
%
%   Input:  D - an n*3 matrix containing the X-, Y- and Z-coordinates of 
%               the points on the original curve.
%           U - an n*1 vector containing the parametrised coordinates of 
%               the points.
%           N - a scalar representing the number of control points
%               required.
%           K - a scalar representing the order of the B-spline curve.
%           V - a 1*m vector containing the knots of the B-spline curve.
%
%   Output: C - an N*3 matrix containing X-, Y- and Z-coordinates of the 
%               control points.
%
% -------------------------------------------------------------------------
% See also: bspbasis, bspchordlparam, bspgenpoints, bspknot, bspcurvefit

U = bspbasis(u, ctrlPNum, K, knotV);

CX=U\data(:,1);
CY=U\data(:,2);
CZ=U\data(:,3);

ctrlp=[CX,CY,CZ];