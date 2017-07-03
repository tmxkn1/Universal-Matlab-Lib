function [R] = FuncICP_rotMat(theta,a,r)

if theta == 0
    R=eye(3);
    return;
end

if r=='d'
    theta=theta*pi/180;
else
end


if a=='x'
    R=[
        1       0                   0;
        0       cos(theta)        -sin(theta);
        0       sin(theta(1))        cos(theta)];
elseif a=='y'
    R=[
        cos(theta)    0       sin(theta);
        0               1       0;
        -sin(theta)   0       cos(theta)];
elseif a=='z'
    R=[
        cos(theta)        -sin(theta)   0 ;
        sin(theta)        cos(theta)    0 ;
        0                   0               1];
else
    
    R=eye(3);
    
    
end

