function [ trans,d_avg ] = FuncICP_Main_ori(P_mp,P_surf,tou,d_min,k_max,n_o_app)
%ICP_func finds the transformation matrix for 
%   Detailed explanation goes here

%% The ICP startes here, where some initial parameters are defined
[~,N]=size(P_mp);


% Center of gravity of both the surface and the measured points are
% calculated and the vector between them
mu_p_mp=FuncICP_centreOfG(P_mp);
mu_p_surf=FuncICP_centreOfG(P_surf);

% The length scale is found
l_scale=zeros(1,3);
for i=1:3
l_scale(i)=max(P_surf(i,:))-min(P_surf(i,:));
end
L_scale=max(l_scale);

e=1/0;
k_count=1;
k_to_t=1;
converge=0;
k_next=-1;
k=0;
d=[1/0];
while d(end)>=d_min && k_count<=n_o_app
    tic
 d=[1/0];             % d is reset for every aproach
 Rot_init(:,:,k_count)=eye(3);  % Initial rotation matrix created
 P=P_mp;            % Measured points are reseat for every aproach
 
 
    if k_count==1    % Initial alignment error is not changed if this is the first approach
    T_rand(:,k_count)=[0;0;0];
    else
 
        if k_count>=50 && k_count<=(50+47)
            P=rot_mat(pi,'x','r')*P;
            Rot_init(:,:,k_count)=Rot_init(:,:,k_count)*rot_mat(pi,'x','r');
        elseif k_count>=98 && k_count<=(98+47)
            P=rot_mat(pi,'y','r')*P; 
            Rot_init(:,:,k_count)=Rot_init(:,:,k_count)*rot_mat(pi,'x','r');
        elseif k_count>(98+47)
            P=rot_mat(pi,'z','r')*P;
            Rot_init(:,:,k_count)=Rot_init(:,:,k_count)*rot_mat(pi,'x','r');
        else
     
        end 

        theta_x=pi/4*((k_count-1)/3+2/3);
        R_x=FuncICP_rotMat(theta_x,'x','r');

        theta_y=pi/4*((k_count-2)/3+2/3);
        R_y=FuncICP_rotMat(theta_y,'y','r');

        theta_z=pi/4*((k_count-3)/3+2/3);
        R_z=FuncICP_rotMat(theta_z,'z','r');

        T_rand(:,k_count)=[(rand-0.5)*L_scale;(rand-0.5)*L_scale;(rand-0.5)*L_scale];
        
        P=P-[mu_p_mp(1)*ones(1,N);mu_p_mp(2)*ones(1,N);mu_p_mp(3)*ones(1,N)];
        
        if k_to_t==1
            P=R_x*P;
            Rot_init(:,:,k_count)=Rot_init(:,:,k_count)*R_x;
            k_to_t=k_to_t+1;
        elseif k_to_t==2
            P=R_y*P;
            Rot_init(:,:,k_count)=Rot_init(:,:,k_count)*R_y;
            k_to_t=k_to_t+1;
        else
            P=R_z*P;
            Rot_init(:,:,k_count)=Rot_init(:,:,k_count)*R_z;
            k_to_t=1;
        end
    
        P=P+[(T_rand(1,k_count)+mu_p_surf(1))*ones(1,N);(T_rand(2,k_count)+mu_p_surf(2))*ones(1,N);(T_rand(3,k_count)+mu_p_surf(3))*ones(1,N)];
 
    end
    
while converge<2 && d(end)>d_min && k<k_max
    k=k+1;
   
    
    % The following function finds the corresponding closest points   
    [ X,~,~] = FuncICP_closestPoint(P,P_surf);
    
   
    % The center of grvity is found for both sets of points
    mu_p=FuncICP_centreOfG(P);
    mu_x=FuncICP_centreOfG(X);
    
    % The cross covariance matrix is calculated   
    [Sigma]=FuncICP_crossCovar(P,X,mu_p,mu_x);
    % The anti symmetric part of the cross covariance matrix is found    
    Theta=Sigma-Sigma';
    % The cyclic components are found from the anti symmetric matrix
    Gamma=[Theta(2,3) Theta(3,1) Theta(1,2)]';
    % The Q matrix is constructed   
    B=Sigma+Sigma'-trace(Sigma)*eye(3);
    Q=[trace(Sigma) Gamma';Gamma B];
    
    % Eigen vectors and eigenvalues are found
    [eig_Vec,eig_Val]=eig(Q);
    
    % The maximum eigenvalue is found
    C_1=max(eig_Val);
    [~,C_2]=max(C_1);
    % The unit quaternion vector is defined as the eigenvector
    % correspoinding to the maximum eigenvalue
    q_R=eig_Vec(:,C_2);
    q=q_R;
    % The optimum rotationa matrix is  found
    R(:,:,k)=[q(1)^2+q(2)^2-q(3)^2-q(4)^2 2*(q(2)*q(3)-q(1)*q(4)) 2*(q(2)*q(4)+q(1)*q(3));
       2*(q(2)*q(3)+q(1)*q(4)) q(1)^2-q(2)^2+q(3)^2-q(4)^2 2*(q(3)*q(4)-q(1)*q(2));
       2*(q(2)*q(4)-q(1)*q(3)) 2*(q(3)*q(4)+q(1)*q(2)) q(1)^2-q(2)^2-q(3)^2+q(4)^2];
   
   % The translation vector is found
   q_T(:,k)=mu_x-R(:,:,k)*mu_p;
   
   % The measured points are rotated
    P=R(:,:,k)*P;
   
   % The measured points are translated
    Q_T=[ones(1,N)*q_T(1,k);ones(1,N)*q_T(2,k);ones(1,N)*q_T(3,k)];
    P=P+Q_T;
    
  
    d(k)=sum(sum((P-X).^2).^0.5)/N;  % The mean error is found

    % Convergence is checked
    if k==1
      e=d(k);
    else
      e=abs(d(k)-d(k-1))  ;           
    end
    k_conv=k;
    if e<tou
       converge=converge+1 ;
       if k_conv==(k_next+1) && converge==2
           
       elseif converge==2
           converge=0;
           k_next=k_conv;           
       else
           k_next=k_conv;
       end
       
    end

end

k_col(k_count)=k;       % Number of iterations of the aproach collected
% % Transolation matrix of the aproach collected
Trans_collect{k_count}=q_T(:,:);
% % Rotational arrey of the aproach collected
Rotation_collect{k_count}=R(:,:,:);
% The final mean error of the aporach collected
d_col(k_count)=d(end);
% Parameters for next aproach set
k_count=k_count+1
k=0;
e=1/0;
converge=0;
toc
end

%% The following section calcualtes the best transformation matrix found
%Best aproach is found by the mean error
[d_avg, I]=min(d_col);
% Optimal transforamtion matrix initiated
trans(:,:)=eye(4);

num_of_iter=k_col(I);
Trans=cell2mat(Trans_collect(I));
Rotation=cell2mat(Rotation_collect(I));
for i=0:(num_of_iter-1)
    
    trans(:,:)=trans(:,:)*[eye(3) Trans(:,(num_of_iter-i));0 0 0 1];
    
    trans(:,:)=trans(:,:)*[Rotation(:,:,(num_of_iter-i)) zeros(3,1);0 0 0 1];
    
end

if I==1
    
else
    trans(:,:)=trans(:,:)*[eye(3) (T_rand(:,I)+ mu_p_surf(1:3,1));0 0 0 1];
    trans(:,:)=trans(:,:)*[Rot_init(:,:,I) zeros(3,1);0 0 0 1];
    trans(:,:)=trans(:,:)*[eye(3) mu_p_mp(1:3,1)*(-1);0 0 0 1];
end


end

