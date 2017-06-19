function a = angle2d(varargin)
%ANGLE2D angle between two 2-D vectors in (-pi, pi] range.
%
% A = angle2d(V,W) returns angle A in radians between vectors V and W.
%
% A = angle2d(...,'Degree') returns angle A in degrees.
%
% by Zhengyi Jiang, University of Manchester, 2017

narg = numel(varargin);
isdeg = 0;

if narg<2
    error('Not enough input arguments.');
end

if narg>3
    error('Too many input arguments.');
end

if narg == 3
    s = varargin{3};
    if strcmp(s,'Degree')
        isdeg = 1;
    else
        error('Invalid parameter ''%s''. A valid parameter is ''Degree''.',s)
    end
end
    
v1 = varargin{1};
v2 = varargin{2};

v1 = v1/norm(v1);
v2 = v2/norm(v2);

a = atan2(v1(1)*v2(2) - v1(2)*v2(1), dot(v1,v2));

if isnan(a)
    a = 0;
end

if isdeg
    a = a*180/pi;
end