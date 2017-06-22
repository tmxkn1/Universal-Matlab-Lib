function R = ptaxisrot (varargin)
%PTAXISROT calculates transformation matrices for rotations about a given 
% axis by a specified angle
%
% R = ptaxisrot (V, P, A) calculates transformation matrices R for
% rotations about an axis defined by a direction vector V and a point P by
% angle A.
%
% required file: VECROTMAT
%
% by Zhengyi Jiang, University of Manchester, 2017

if nargin < 3
    error('Not enough input arguments.');
end

if nargin > 4
    error('Too many input arguments.');
end

v = varargin{1}; p = varargin{2}; a = varargin{3};

if nargin == 4
    if strcmp(varargin{4},'Degree')
        a = a*pi/180;
    else
        error('Invalid parameter ''%s''. The only valid parameter is ''Degree''.', s);
    end
end

% normalise the direction vector in case it hasn't been done yet
v = v(:)/norm(v);
p = p(:);

% transform axis to x-axis
r = vecrotmat(p, v+p, [0 0 0], [1 0 0]);

% rotation matrix about x
rx = [1 0 0 0; 0 cos(a) sin(a) 0; 0 -sin(a) cos(a) 0; 0 0 0 1];

R = r * rx * r^-1; 