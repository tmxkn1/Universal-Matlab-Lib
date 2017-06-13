function varargout = vecrotmat(a1, a2, b1, b2)
%VECROTMAT calculates a matrix that transforms one segment to the other 
% such that the first segment overlaps the second and the two start points
% are also overlapping.
%
% R = vecrotmat(a1, a2, b1, b2) calculates a matrix, R, that transforms 
% segment A1-A2 to segment B1-B2.
% 
% [R, T1, T2, R1, R2, R3] = vectotmat(...) also returns all individual
% translations and rotations.
%
% Running the function without input for a demo, i.e. vecrotmat.
%
% Example use of the result:
% [ (b1 - b2) / norm(b1 - b2), 1 ] = [ (a1 - a2) / norm(a1 - a2), 1 ] * R
%
% Reference:
% angle between two vectors:
% https://uk.mathworks.com/matlabcentral/answers/16243-angle-between-two-vectors-in-3d
% rotation matrix (orthogonal):
% https://en.wikipedia.org/wiki/Rotation_matrix

if nargin == 0
    demo = 1;
else
    demo = 0;
end

a1 = a1(:)';
a2 = a2(:)';
b1 = b1(:)';
b2 = b2(:)';

if demo
    a1 = sign((randi(2,1,3)-1.5)).*randi(10,1,3);
    a2 = sign((randi(2,1,3)-1.5)).*randi(10,1,3);
    b1 = sign((randi(2,1,3)-1.5)).*randi(10,1,3);
    b2 = sign((randi(2,1,3)-1.5)).*randi(10,1,3);
end

% In this example lineA is going to be transformed to where lineB is, so they 
% overlap.
% lineA is defined by two points a1 and a2 and
% lineB is defined by two points b1 and b2.

%% Find the direction vectors of the two lines, va and vb
va = a1-a2;
vb = b1-b2;

%% Translations
% Find the translations Ta and Tb which move a1 and b1 to the origin
% respectively
Ta = -a1;
Tb = -b1;


% Find the orthogonal translations T1 and T2
T1 = [[1 0 0; 0 1 0; 0 0 1; Ta] [0; 0; 0; 1]]; % a1 * T1 = a1 + Ta
T2 = [[1 0 0; 0 1 0; 0 0 1; Tb] [0; 0; 0; 1]]; % b1 * T2 = b1 + Tb

%% Now we have 
% a1 * T1 = b1 * T2, 
% i.e. 
% b1 = a1 * T1 * T2^-1.

% By applying the transformation T1 * T2^-1 to lineA, we have lineA', 
% which intersects lineB at b1.

% Although it is much simpler to define the translation as T = b1-a1 and 
% by having lineA + T we can achieve the same result , there is a good
% reason why we want to use T1 and T2, i.e. the translations that moves 
% lineA and lineB to the origin.

% Let's assume we have lineA' = lineA + T. Next, we need to rotate lineA'
% about b1 to make it overlap lineB. This can be done by a single rotation 
% about the axis of symmetry. However, the solution is difficult to find, 
% because we only have the equations for rotating about the X-, Y- and 
% Z-axis.

% A better method is to move both lines to the orgin:
% lineA' = lineA * T1
% lineB' = lineB * T2
% then rotate them about the three axes. This is why we have to calculate 
% T1 and T2.

% To calculate the rotations, we need to use direction vectors, which have
% been already calculated as va and vb.

% The transformation can be done with two rotations:
% Let us assume the two vectors, v1 and v2, are overlapping and we want to 
% seperate them by rotating v1 about the axes.

% STEP 1, we rotate v1 about the Z-axis, and we have v1' and v2. Since the
% rotation is about the Z-axis, the angle between any of the vectors and the
% Z-axis is not changed.

% STEP 2, we rotate v1' about the Y-axis, so it becomes v1''.

% Now let's assume v1'' is va and v2 is vb. We want to rotate va so it
% overlaps vb. We reverse the above steps.

% STEP 1, rotate va about the Y-axis and we have va'. This rotation must
% achieve a result where the angle between va' and the Z-axis equals the 
% angle between vb and the Z-axis.

% STEP 2, rotate va about the Z-axis until it overlalps vb.

% This is only a simple example. In real world problems, we can't ensure at
% the end of STEP 1, i.e. one rotation about the Y-axis can achieve what we
% need at the start of STEP 2. 
% For example: 
% va = [0.98, 0.13, 0.13] and vb = [0.001, 0.001, 0.999].

% Normally, we need two rotations at STEP 1.
% The first rotation is to rotate va about the X- or Z-axis until it lies  
% on the XZ-plane: 
% va' ~ [0.9989, 0, 0.0011]. 
% Next, va' can be rotated about the Y-axis to achieve the position 
% required: 
% va'' ~ [0.0015, 0, 0.999]


% Having these knowledge, we can now move on.

%% STEP 1.1
% Rotate va to the XZ-plane.
% Find the axial angle between the vector and XZ-plane, i.e. the angle
% between the project of va and the x-axis.

% Find the projection of va on the XY-plane, reduced to 2d.
pa = va(1:2);

% Find the angle between the two projection.
ra = angle2d(pa,[1,0]);

% Find the first rotation matrix R1.
R1 = [cos(ra) sin(ra) 0; -sin(ra) cos(ra) 0; 0 0 1];

% Transform va.
va = va * R1;

% Make R1 orthogonal.
R1 = [[R1; [0 0 0]] [0; 0; 0; 1]];

if demo
    figure;hold on; axis equal
    plot3([0 va(1)],[0 va(2)],[0 va(3)],'-o')
    plot3([0 vb(1)],[0 vb(2)],[0 vb(3)],'-^')
    title('R1');legend('va', 'vb');grid on;axis tight;
    xlabel x; ylabel y; zlabel z
end

%% STEP 1.2
% Find the angles between the Z-axis and both vectors, reduced to 2d. 
% Rotate va about the Y-axis by the difference between the angles.
angle1 = atan2(norm(cross([0 0 1],va)),dot([0 0 1],va));
angle2 = atan2(norm(cross([0 0 1],vb)),dot([0 0 1],vb));
ra = angle2 - angle1;

% Find the second rotation matrix R2.
R2 = [cos(ra) 0 -sin(ra); 0 1 0; sin(ra) 0 cos(ra)];

% Transform va.
va = va * R2;

% Make R2 orthogonal.
R2 = [[R2; [0 0 0]] [0; 0; 0; 1]];

if demo
    figure;hold on; axis equal;
    plot3([0 va(1)],[0 va(2)],[0 va(3)],'-o')
    plot3([0 vb(1)],[0 vb(2)],[0 vb(3)],'-^')
    title('R2');legend('va', 'vb');grid on;axis tight;
    xlabel x; ylabel y; zlabel z
end

%% Find the third rotation
% Rotate va about the Z-axis, so that va overlaps vb.
% Find the XY-plane projections of both vectors, reduced to 2d.
pa = [va(1) va(2)];
pb = [vb(1) vb(2)];

% Find the angle between the two projections.
ra = angle2d(pa,pb);

% Find the third rotation matrix R3.
R3 = [cos(ra) sin(ra) 0; -sin(ra) cos(ra) 0; 0 0 1];

% Transform va.
va = va * R3;

% Make R3 orthogonal.
R3 = [[R3; [0 0 0]] [0; 0; 0; 1]];

if demo
    figure;hold on; axis equal;
    plot3([0 va(1)],[0 va(2)],[0 va(3)],'-o')
    plot3([0 vb(1)],[0 vb(2)],[0 vb(3)],'-^')
    title('R3');legend('va', 'vb');grid on;axis tight;
    xlabel x; ylabel y; zlabel z
end

%% Combine the translation and rotations
% Combine rotations and translations.
R = T1*R1*R2*R3*T2^-1;
varargout = {R, T1, T2, R1, R2, R3};

%% Check your results
if demo
    a1s = [ a1, 1 ];
    a2s = [ a2, 1 ];
    a1s = a1s * R;
    a2s = a2s * R;
    
    figure;hold on; axis equal;
    plot3([a1(1) a2(1)],[a1(2) a2(2)],[a1(3) a2(3)])
    plot3([b1(1) b2(1)],[b1(2) b2(2)],[b1(3) b2(3)],'-^')
    plot3([a1s(1) a2s(1)],[a1s(2) a2s(2)],[a1s(3) a2s(3)],'-o')
    title('T1*R1*R2*R3*T2^-1');legend('a1a2', 'b1b2', 'a1a2 transformed');grid on;axis tight;
    xlabel x; ylabel y; zlabel z
    varargout = [varargout, [a1;a2], [b1;b2]];
end

function a = angle2d(v1,v2)
a = atan2(v1(1)*v2(2) - v1(2)*v2(1), dot(v1,v2));