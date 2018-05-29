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
    disp('line line has not been implemented');
elseif isequal(types, [2,3])
    % line plane
    [obj, condi] = planeline(inputobj{1}, inputobj{2});
    disp('line plane')
elseif isequal(type, [3,3])
    % plane plane
    disp('plane plane')
elseif iseuqla(type, [2,4])
    % line circle
    disp('line circle')
end


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

I=nan(1, line.Dimension);
u = line.DirectionVector;
w = line.Point.Value - plane.Point.Value;
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
I = line.Point.Value+ sI.*u;
C = 1;













