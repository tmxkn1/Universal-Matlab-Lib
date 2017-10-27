function [c, n, r, e, pfit] = circfit3d(p)
% CIRCFIT3D fits a best circle to a set of points in 3D space.
%
% [C, N, R] = circfit3d(P) returns the centre C, normal N, and radius R, of
% the best circle fitted to a set of points P(x,y,z). P is an Nx3 matrix.
%
% [..., E] = circfit3d(P) also returns fitting error.
%
% [..., PFIT] = circfit3d(P) also returns points calculated on the fitted
% circle using P(:,1).

pmean = mean(p);
np = length(p);
pcentred = p - ones(np, 1)*pmean;

% find the best plane
[~, ~, n]= svd(pcentred);
n = n(1:3);
% rotate and project the points
R = rot(n);
prot = [pcentred ones(np,1)]*R;
prot(:,4) = [];
s = sum(abs(prot),1);
[~,I] = min(s);
d1 = abs(5-I*2);
d2 = 6-d1-I;
% best circle
[xc, yc, r, ~, e] = circfit(prot(:, d1), prot(:, d2));
c = [xc, yc, 0] + pmean;

if nargout == 5
    x = linspace(min(prot(:,d1)),max(prot(:,d1)),np);
    y = sqrt( r^2 - (x-xc).^2 ) + yc;
    y2 = -sqrt( r^2 - (x-xc).^2 ) + yc;
    
    if mean( abs(y-prot(:,d2)) ) > mean( abs(y2-prot(:,d2)) )
        y = y2;
    end
    prot(:,d1) = x';
    prot(:,d2) = y';
    pfit = [prot ones(np,1)] * R^-1;
    pfit = pfit(:,1:3) + pmean;    
end


function R = rot(v)
% project v on the xy plane and find the angle between the projection and
% the x axis
a = atan2(norm(cross([v(1),v(2),0],[1,0,0])), dot([v(1),v(2),0],[1,0,0]));
% rotate v to x-z plane
R1 = rotz(a)^-1;
v = v*R1;
% find the angle between rotated v and the z-axis
a = atan2(norm(cross(v,[0,0,1])),dot(v,[0,0,1]));
% rotation matrix
R2 = roty(a)^-1;
R1 = [[R1; [0 0 0]] [0; 0; 0; 1]];
R2 = [[R2; [0 0 0]] [0; 0; 0; 1]];
R = R1*R2;