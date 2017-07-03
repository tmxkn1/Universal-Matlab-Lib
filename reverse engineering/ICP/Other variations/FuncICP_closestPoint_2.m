function [ closest_points,distance,places] = FuncICP_closestPoint_2(meas_points,surf_points)
%   This function calculates the closest points at the surface f_surf by
%   differanciating


[M_m,N_m]=size(meas_points);                     % Establish how many point there are
[M_s,N_s]=size(surf_points);

x_surf=surf_points(1,:);
y_surf=surf_points(2,:);

if M_s==2
    z_surf=zeros(1,N_s);
else
    z_surf=surf_points(3,:);
end
x_meas=meas_points(1,:);
y_meas=meas_points(2,:);

if M_m==2
    z_meas=zeros(1,N_m);
else
    z_meas=meas_points(3,:);
end
closest_points=zeros(M_m,N_m);
distance=zeros(N_m,1);
places=zeros(N_m,1);
surfP = [x_surf;y_surf;z_surf];
%b_temp_ones = ones(1,N_s);
measP = [x_meas;y_meas;z_meas];
parfor i=1:N_m
    %b_temp = measP(:,i)*b_temp_ones;
    
    dist=sqrt(bsxfun(@minus, x_surf, measP(1,i)).^2+bsxfun(@minus, y_surf, measP(2,i)).^2+bsxfun(@minus, z_surf, measP(3,i)).^2);
    if y_meas(i)<-26400
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

