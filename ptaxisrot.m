function R = ptaxisrot (v, p, a)
%PTAXISROT calculates transformation matrices for rotations about a given 
% axis by a specified angle
%
% R = ptaxisrot (V, P, A) calculates transformation matrices R for
% rotations about an axis defined by a direction vector V and a point P by
% angle A.
%
% requires VECROTMAT

% normalise the direction vector in case it hasn't been done yet
v = v(:)/norm(v);
p = p(:);

% transform axis to x-axis
r = vecrotmat(p, v+p, [0 0 0], [1 0 0]);

% rotation matrix about x
rx = [1 0 0 0; 0 cosd(a) sind(a) 0; 0 -sind(a) cosd(a) 0; 0 0 0 1];

R = r * rx * r^-1; 