function [obj, condi] = geointersect ( varargin )

if nargin < 2
    error('Not enough input arguments.')
end

types = [int32(varargin{1}.Type) int32(varargin{2}.Type)];
inputobj = {varargin{1} varargin{2}};
[types, id] = sort(types);
inputobj = inputobj(id);

if isequal(types, [2,2])
    % line line
    if inputobj{1}.Dimension == 3
        [obj, condi] = lineline3D(inputobj{1}, inputobj{2});
    elseif inputobj{1}.Dimension == 2
        [obj, condi] = lineline2D(inputobj{1}, inputobj{2});
    end
        
elseif isequal(types, [2,3])
    % line plane
    [obj, condi] = planeline(inputobj{1}, inputobj{2});
elseif isequal(types, [3,3])
    % plane plane
    disp('plane plane')
elseif isequal(types, [2,4])
    % line circle
    [obj, condi] = linecirc(inputobj{1}, inputobj{2});
end

function [I, C]=lineline3D(l1, l2)
%oUTPUTS:
%
% I is the point of intersection 
%
% C indicates intersection condition: 
%   0 => no intersection
%   1 => parallel or colinear
%   2 => intersects

I = nan(1, 3);
C = 0;
g = l1.Point - l2.Point;
h = norm(cross(l1.DirectionVector, g));
if h == 0
    return % no intersection
end
k = norm(cross(l1.DirectionVector, l2.DirectionVector));
if k == 0
    C = 1; % parallel
    return
end
C = 2; % intersection found
t = h/k * l2.DirectionVector;
I1 = l2.Point + t;
if l1.isOnLine(I1)
    I = I1;
    return
end
I = l2.Point - t;    

function [I, C]=lineline2D(l1, l2)
%oUTPUTS:
%
% I is the point of intersection 
%
% C indicates intersection condition: 
%   1 => parallel or colinear
%   2 => intersects

t = (l2.Point - l1.Point) / [l1.DirectionVector; -l2.DirectionVector];

if any(abs(t) == [inf inf])
    C = 1;
    I = nan(1, 2);
    return;
end

C = 2;
I = t(1) * l1.DirectionVector + l1.Point;

function [I, C]=planeline(line, plane)
%Outputs: 
%
% I is the point of intersection 
%
% C indicates intersection condition: 
%   0 => parallel (no intersection) 
%   1 => the plane intersects the line at the unique point I 
%   2 => the line lies in the plane 
%
% This function is adopted from "plane_line_intersection" written by Nassim
% Khaled, Wayne State University, Research Assistant and Phd candidate

I = nan(1, line.Dimension);
u = line.DirectionVector;
w = line.Point - plane.Point;
D = dot(plane.NormalVector, u);
N = -dot(plane.NormalVector, w);

if abs(D) < realmin('single') % The segment is parallel to plane
    if N == 0 % The segment lies in plane
        C=2;
        return
    else
        C=0; %no intersection
        return
    end
end

%compute the intersection parameter
sI = N / D;
I = line.Point+ sI.*u;
C = 1;

function [I,C] = linecirc(line,circ)
%Outputs: 
%
% I is a 2-by-n matrix, each row representing one point of intersection 
%
% C indicates intersection condition: 
%   0 => no intersection 
%   1 => the line intersects the circle at one point
%   2 => the line intersects the circle at two points
%
% by Zhengyi Jiang, ELE Advanced Technologies Ltd, 2017

I = [];
C = 0;

% Solve t^2 * (dot(Ld,Ld)) + 2t*(dot(-C,Ld)) + (dot(C,C) - r^2 ) = 0 for t,
% where Ld is the direction vector of the line, C is the centre of the
% circle, r is the radius, and t = PoI(x,y,z...)/Ld(x,y,z...)
d = line.Point-circ.Centre;
a = dot(line.DirectionVector, line.DirectionVector);
b = 2*dot(d,line.DirectionVector);
c = dot(d, d) - circ.Radius^2;
t = b^2-4*a*c;

% no intersection
if t < 0
    return
end

t = sqrt(t);
t1 = (-b - t)/2/a;
t2 = (-b + t)/2/a;

if t1==t2
    C = 1;
else
    C = 2;
end

I = [t1*line.DirectionVector + line.Point
    t2*line.DirectionVector + line.Point];