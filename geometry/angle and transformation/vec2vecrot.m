function R = vec2vecrot(v1,v2)
%VEC2VECROT calculates a transformation matrix that rotates and translates
% vector 1 so that vector 1 is colinear with vector 2.
% 
% Input:
%   v1, v2      both must either be 1-by-3 vectors or 2-by-3 matrices:
%
%   if 1-by-3   v1 and v2 define the vectors.
%   if 2-by-3   each row in v1 and v2 defines a point, the two points in v1
%               and v2 define the vectors respectively.
% Output:
%   R           if inputs are 1-by-3, R is a 3-by-3 transformation matrix;
%               otherwise, it is a 4-by-4 homogeneous transformation
%               matrix.
%               Geometric data to be transformed must be multiplied by R, 
%               i.e. Data*R.      

if size(v1,1)==2
    T1 = eye(4);
    T2 = T1;
    T1(4,1:3) = - v1(1,:);
    T2(4,1:3) = v2(1,:);
    R = vec2vecrot(v1(2,:)-v1(1,:), v2(2,:)-v2(1,:));
    R = T1*tmat2hom(R^-1)*T2;
    return
end

% normalise vectors
v1 = v1/norm(v1,'fro');
v2 = v2/norm(v2,'fro');

% axis of rotation
ax = cross(v1,v2);

% Rodrigues rotation
R = rodriguesrot(ax, norm(ax,'fro'), dot(v1,v2))^-1;