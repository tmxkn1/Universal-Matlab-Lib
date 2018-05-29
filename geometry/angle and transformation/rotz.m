function r = rotz(varargin)
%ROTX rotation matrix about z-axis
% 
% R = rotx(A) returns R representing the rotation about the z-axis by angle
% A in radians.
%
% R = rotx(A,'Degree') returns R representing the rotation about the z-axis
% by angle A in degrees.
%
% R = rotx(...,'Homo') returns R in the homogeneous form.
%
% R is a rotation matrix that transforms P1 via a rotation to P2: P2 = P1*R
%
% by Zhengyi Jiang, ELE Advanced Techonologies Ltd., 2017

isHomo = 0;

if nargin < 1
    error('Not enough input arguments.');
end

if nargin > 4
    error('Too many input arguments.');
end

a = varargin{1};

if nargin > 1
    for i = 2:nargin
        s = varargin{i};
        if ~ischar(s)
            error('Parameter %d is invalid. Valid parameters are ''Degree'' and ''Homo''.', i)
        end
        if strcmp(s, 'Degree')
            a = a/180*pi;
        elseif strcmp(s, 'Homo')
            isHomo = 1;
        else
            error('Invalid parameter ''%s''. Valid parameters are ''Degree'' and ''Homo''.', s);
        end
    end
end

r = [cos(a) sin(a) 0; -sin(a) cos(a) 0; 0 0 1];

if isHomo
    r = [ [ r ; 0 0 0] [ 0; 0; 0; 1 ] ];
end