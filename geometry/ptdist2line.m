function d=ptdist2line(p, v, c)

% vector cp
cp = c-p;
% area of the parallelogram with sides given by p-c and v
area = norm(cross(cp, v),'fro');
% height of the parallelogram
d = area/norm(v,'fro');