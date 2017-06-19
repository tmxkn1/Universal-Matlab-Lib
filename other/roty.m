function r = roty(varargin)
%ROTX rotation matrix about x-axis
% 
% R = rotx(A) returns R representing the rotation about the y-axis by angle
% A in radians.
%
% R = rotx(A,'Degree') returns R representing the rotation about the y-axis
% by angle A in degrees.
%
% R = rotx(...,'Homo') returns R in the homogeneous form.
%
% R is a rotation matrix that transforms P1 via a rotation to P2: P2 = P1*R
%
% by Zhengyi Jiang, ELE Advanced Techonologies Ltd., 2017
narg = numel(varargin);

isDeg = 0;
isHomo = 0;

if narg < 1
    error('Not enough input arguments.');
end

if narg > 4
    error('Too many input arguments.');
end

a = varargin{1};

if narg > 1
    for i = 2:narg
        s = varargin{i};
        if ~ischar(s)
            error('Parameter %d is invalid. Valid parameters are ''Degree'' and ''Homo''.', i)
        end
        if strcmp(s, 'Degree')
            isDeg = 1;
            continue;
        elseif strcmp(s, 'Homo')
            isHomo = 1;
        else
            error('Invalid parameter ''%s''. Valid parameters are ''Degree'' and ''Homo''.', s);
        end
    end
end


if isDeg
    a = a/180*pi;
end

r = [cos(a) 0 -sin(a); 0 1 0; sin(a) 0 cos(a)];

if isHomo
    r = [ [ r ; 0 0 0] [ 0; 0; 0; 1 ] ];
end