function R = rodriguesrot(ax, s, c)
%RODRIGUESROT calculates a rotation by a given angle about a given axis.
%
% R = rodriguesrot(AX, A) calculates a rotation matrix R from the rotary
% axis AX and the angle of rotation A.
%
% NOTE: rotation follows the left hand rule.
%
% R = rodriguesrot(AX, S, C) calculates a rotation matrix R from the rotary
% axis AX and sin, S, and cosin, C, of the angle.
%
% Use:
%
% To rotate point P by angle A, about axis AX:
%           P_new = rodriguesrot(AX,A) * P' 

if nargin == 2
    c = cos(s);
    s = sin(s);
end

% in the case where two vectors are colinear/parallel.
if all(ax==0)
    R = eye(3);
    return
end

% normalise the direction vector in case it hasn't been done yet
ax = ax(:)/norm(ax, 'fro');

% calculate K
K = [0 -ax(3) ax(2); ax(3) 0 -ax(1); -ax(2) ax(1) 0];

% calculate rotation matrix
R = eye(3) + s * K + (1-c) * K^2;