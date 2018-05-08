function p = pt2hom(p)
%PT2HOM converts input coordinates to homogeneous coordinates.
% The input must be a row-wise matrix, each row is a set of coordinate.

p = [p ones(size(p,1),1)];