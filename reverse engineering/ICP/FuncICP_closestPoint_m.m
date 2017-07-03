function [ closest_points,mindist,match] = FuncICP_closestPoint_m(Q, P)
% Q target point set
% P floating point set
% For each point in P, find its closest point in Q.
% i.e. P covers a smaller area.

testresult = license('test', 'Statistics_Toolbox');

if size(Q,1) ~= 2
    Q = Q';
end
if size(P,1) ~= 2
    P = P';
end

if testresult == 1
    kdOBJ = KDTreeSearcher(Q);
    [match, mindist] = knnsearch(kdOBJ,P);
    closest_points = Q(match,:);
    closest_points = closest_points';
else
    [ closest_points,mindist,match] = match_distance(P,Q);
end

function [ closest_points,distance,places] = match_distance(P,Q)
% Q target point set
% P floating point set

[M_m,N_m]=size(P);
[M_s,N_s]=size(Q);

x_surf=Q(1,:);
y_surf=Q(2,:);

if M_s==2
    z_surf=zeros(1,N_s);
else
    z_surf=Q(3,:);
end
x_meas=P(1,:);
y_meas=P(2,:);

if M_m==2
    z_meas=zeros(1,N_m);
else
    z_meas=P(3,:);
end
closest_points=zeros(M_m,N_m);
distance=zeros(N_m,1);
places=zeros(N_m,1);
b_temp_ones = ones(1,N_s);
measP = [x_meas;y_meas;z_meas];
for i=1:N_m
    b_temp = measP(:,i)*b_temp_ones;
    
    dist=((b_temp(1,:)-x_surf).^2+(b_temp(2,:)-y_surf).^2+(b_temp(3,:)-z_surf).^2).^0.5;
    if y_meas(i)<-2669
        fval=0;
        count_zeros=count_zeros+1;
    else
        [fval,I] = min(dist);
        closest_point=[x_surf(I) y_surf(I) z_surf(I)];
        distance(i)=fval;
        closest_points(:,i)=closest_point;
        places(i)=I;
    end
end