function fitp = bspfit(u, ctrlp, K, knotV)
% Input arguments:
%   u:
%       a row vector of n number of u coordinates.
%
%   ctrlp:
%       a n x 3 matrix; each column is the x-, y- and z-coordinates of all
%       control points, can be obtained with:
%       ctrlp = sub_Bspline_ctrlp(data, u, no_of_ctrl, deg_of_curve, 
%                                   knotV, d1, d2)
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
%   fitp:
%       a no_of_u x 3 matrix; each column is the x-, y- and z-coordinates 
%       of all fitted points

U = bspbasis(u, numel(ctrlp(:,1)), K, knotV);

FX=U*ctrlp(:,1);
FY=U*ctrlp(:,2);
FZ=U*ctrlp(:,3);

fitp=[FX,FY,FZ];