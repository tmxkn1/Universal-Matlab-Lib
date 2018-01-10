function wpr = vectowpr(vector, r)
% VECTOWPR converts a direction vector in Cartesian coordinates to WPR
% representation, which is used by certain industrial robots such as FANUC.
%
% Warning: With this solution, r is not calculated. The output of r is
% zero.
% This method converts a vector into W,P,R with a method involving rotation
% matrix.
%
% WPR = VECTOWPR(V) converts the vector V in Cartesian coordinates to WPR
% representation, with R set to zero.
% 
% WPR = VECTOWPR(V,R) allows user to specify an R value.
%
% =========================================================================
% LOGICS:
% Assume we have an arbitrary unit vector, A, and a unit vector N=[0 0 1].
% N needs to be rotated at most three times about x, y and z axis in
% sequence to be aligned with A.
% The angles of rotations are defined as w, p and r respectively and we
% need to calculate them. This can be done by solving rotation matrices.

% Rotation about x is defined by:
% Rx = [ 1 0 0; 0 cos sin; 0 -sin cos ];
% Rotation about y is defined by:
% Ry = [ cos 0 -sin; 0 1 0; sin 0 cos ];
% Rotation about z is defined by:
% Rz = [ cos sin 0; -sin cos 0; 0 0 1 ];
% Therefore, we have A = N * Rx(w) * Ry(p) * Rz(r).
% Let r = 0, we only need to simplify A = N * Rx(w) * Ry(p).
% So we have A = [ cos(w)*sin(p), -sin(w), cos(p)*cos(w) ].
% Finally: A(1) = cos(w)*sin(p) ... Eq.1
%          A(2) = -sin(w) ......... Eq.2
%          A(3) = cos(p)*cos(w) ... Eq.3
% Either for w or p, we have two solutions taking quadrant into
% consideration. Eventually, we end up with four different pairs of
% solution and we only pick one from them. This will be decided by
% converting each pair of solution back to the vector.
%
% If we are getting an r value, for a specific end effect orientation, we
% can always have A/Rz(r) = N * Rx(w) * Ry(p);

isTesting = 0;

if nargin == 1
    r = 0;
end

%% Find all four solutions
% Normalize vector
A_ori = vector;
if size(A_ori,1) == 3
    A_ori = A_ori';
end
A_ori = A_ori/norm(A_ori);

% A/Rz(r)
Rz = [ cosd(r) sind(r) 0; -sind(r) cosd(r) 0; 0 0 1 ];
A = A_ori/Rz;

% Get two solutions for w by solving Eq.2
w1 = asind(-A(2));
% When sin > 0, the second solution w2 in [-pi, pi] is 180 - w1.
% When sin < 0, the second solution w2 in [-pi, pi] is -180 - w1.
if -A(2) > 0
    w2 = 180 - w1;
else
    w2 = -180 - w1;
end

% Get two solutions for p by solving Eq.1/Eq.3
p1 = atand(A(1)/A(3));
% When tan > 0, the second solution p2 in [-pi, pi] is -180 + p1;
% When tan < 0, the second solution p2 in [-pi, pi] is 180 + p1;
if A(1)/A(3) > 0
    p2 = -180 + p1;
else
    p2 = 180 + p1;
end

w = [w1 w2];
p = [p1 p2];

%% Check solutions
% % We want w > p
% sol = [0 0 0 0];
% c = 1;
% for i = 1:2
%     for j = 1:2
%         if abs(w(i)) > abs(p(j))
%             sol(c) = 1;
%         end
%         c = c+1;
%     end
% end
%
% % Calculate the vector with passed combination of w, p
% c = 1;
% for i = 1:2
%     for j = 1:2
%         if sol(c) == 1
%             s = FuncUSToolPath_wprTransformation(w(i), p(j), 0, 0);
%             if mean(abs(s-A)) < 0.0001
%                 sol(c) = 1;
%             else
%                 sol(c) = 0;
%             end
%         end
%         c = c+1;
%     end
% end
%
% % Format results
% wpr = [w1 p1 0
%     w1 p2 0
%     w2 p1 0
%     w2 p2 0];
% wpr = [sol;sol;sol]'.*wpr;
%
% % Delete 0 0 0
% wpr(~any(wpr,2),:)=[];

% Format results [w p r check]
wprc = [w1 p1 r 999
    w1 p2 r 999
    w2 p1 r 999
    w2 p2 r 999];

% Calculate the vector with passed combination of w, p
for i = 1:4
    s = wprtovec(wprc(i,1), wprc(i,2), r);
    if mean(abs(s-A_ori)) < 0.0001
        wprc(i,4) = 1;
    end
end

if isTesting
    A
    wprc
end

% Introduce a unit vector [1 0 0]. This vector is rotated with all four
% solutions. Angles between the resultant vectors and the original vector
% [1 0 0] are checked. The solution results in a minimum absolute angle is
% selected.
V = [ 1 0 0 ];
V = V/Rz;
for i = 1:4
    if wprc(i,4) == 1
       Vnew =  wprtovec(wprc(i,1), wprc(i,2), r, V);
       wprc(i,4) = abs(atan2(norm(cross(V,Vnew)),dot(V,Vnew)));
    end
end

if isTesting
    wprc
end

I = find(wprc(:,4) == min(wprc(:,4)));
wpr = wprc(I,1:3);