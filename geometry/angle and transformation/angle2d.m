function a = angle2d(varargin)
%ANGLE2D angle between two 2-D vectors in (-pi, pi] range.
%
% A = angle2d(V,W) returns angle A in radians from vector V to W. V and W
% do not need to be normalised.
%
% NOTE: Positive is anti-clockwise.
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

s1 = size(v1);
s2 = size(v2);

if ~isequal(s1,s2)
    if isvector(v1)
        v1 = repmat(v1,s2(1),1);
    elseif isvector(v2)
        v2 = repmat(v2,s1(1),1);
    else
        error('Input arguments should be two vectors, two matrices of the same size or one vector and one matrix')
    end
end

v1 = v1./sqrt(sum(v1.^2,2));
v2 = v2./sqrt(sum(v2.^2,2));

a = atan2(v1(:,1).*v2(:,2) - v1(:,2).*v2(:,1), dot(v1',v2')');

a(isnan(a)) = 0;

if isdeg
    a = a*180/pi;
end