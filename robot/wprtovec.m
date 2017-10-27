function [n_new] = wprtovec(w, p, r, vector)
% WPRTOVEC convert robot w-p-r coordinates to a direction vector
%
% [V] = WPRTOVEC(W, P, R) calculates a direction vector V in
% Cartesian coordinates with given W, P and R coordinates as well as the
% initial vector VI=[0,0,1], which is used by a standard FANUC robot.
%
% [V] = WPRTOVEC(..., VI) allows user to specify VI.
%
% =========================================================================
% VI represents the direction vector in Cartesian coordinates at w=p=r=0.
% =========================================================================

if nargin == 3
    n = [0, 0, 1, 0];
else
    n = [vector, 0];
end

if numel(w)~=numel(p) || numel(w)~=numel(r)
    error('Size mismatch. The sizes of arrays W, P and R must equal.');
end

for i = 1:numel(w)
    
    Rx = [1 0 0 0
        0 cosd(w(i)) sind(w(i)) 0
        0 -sind(w(i)) cosd(w(i)) 0
        0 0 0 1];
    
    Ry = [cosd(p(i)) 0 -sind(p(i)) 0
        0 1 0 0
        sind(p(i)) 0 cosd(p(i)) 0
        0 0 0 1];
    
    Rz = [cosd(r(i)) sind(r(i)) 0 0
        -sind(r(i)) cosd(r(i)) 0 0
        0 0 1 0
        0 0 0 1];
    
    n_new(i,:) = n*Rx*Ry*Rz;
    
end

n_new(:,4) = [];