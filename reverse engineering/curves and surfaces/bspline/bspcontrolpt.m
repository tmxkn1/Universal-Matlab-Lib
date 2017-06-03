function ctrlp = bspcontrolpt(data, u, ctrlPNum, K, knotV)
% Input arguments:
%   data:
%       a n x 3 matrix; each column is the x-, y- and z-coordinates of
%       data points.
%
%   u:
%       a row vector of n number of u coordinates.
%
%   ctrlPNum:
%       a scalar, equals the number of control points.
%
%   curveDeg:
%       a scalar, equals the degree of the curve.
%
%   knotV:
%       a row vector of no_of_ctrl number of knot vectors, can be obtained
%       by [knotV, d1, d2] = sub_Bspline_knotv(no_of_ctrl, deg_of_curve).
%
%   d1 and d2:
%       matrices of denominators of basis functions, can be obtained along
%       with knotV.
%
% Output arguments:
%   ctrlp:
%       a no_of_ctrl x 3 matrix; each column is the x-, y- and
%       z-coordinates of all control points.

U = bspbasis(u, ctrlPNum, K, knotV);

CX=U\data(:,1);
CY=U\data(:,2);
CZ=U\data(:,3);

ctrlp=[CX,CY,CZ];