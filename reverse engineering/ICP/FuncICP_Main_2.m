function [ trans,d_avg, icp ] = FuncICP_Main_2(P_mp,P_surf,varargin)
% P_mp smaller set of points
inp = inputParser;
inp.addRequired('P_mp', @(x)isreal(x) && size(x,1) == 3);
inp.addRequired('P_surf', @(x)isreal(x) && size(x,1) == 3);

inp.addOptional('k_max', 100, @(x)x > 0);

inp.addParameter('InitialTrans', eye(4), @(x)isreal(x) && size(x) == [4,4]);
inp.addParameter('ConvergeCrit', 0.000001, @(x)x >= 0);
inp.addParameter('MinE', 0.005, @(x)x >= 0);
inp.addParameter('NoApproach', 1, @(x)x > 0);

validMinimize = {'SVD','Quaternion'};
inp.addParameter('Minimize', 'Quaternion', @(x)any(strcmpi(x,validMinimize)));

inp.parse(P_mp,P_surf,varargin{:});
arg = inp.Results;
clear('inp');


%%
terminationInfo = {'Minimum error reached!','Convergence reached!',...
    'Maximum number of iteration...' };
rot_map = [180,0,0; 0,180,0; 0,0,180;];
%     90,0,0; 0,90,0; 0,0,90; 90,90,0; 90,0,90; 0,90,90; 90,90,90;
%     -90,0,0 ; 0,-90,0 ; 0,0,-90 ; -90,-90,0 ; -90,0,-90; 0,-90,-90; -90,-90,-90;];

if arg.NoApproach > size(rot_map,1)
    arg.NoApproach = size(rot_map,1)+1;
end

%% The ICP startes here, where some initial parameters are defined
[~,N]=size(P_mp);

if arg.InitialTrans==0
    arg.InitialTrans=eye(4);
end

% Center of gravity of both the surface and the measured points are
% calculated and the vector between them
mu_p_mp=P_mp;
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
d=1/0;
while d(end)>=arg.MinE && k_count<=arg.NoApproach
    
    clear qr qt
    timeCost = 0;
    terminateCondi_ = 0;
    nextProgressReport = 0;
    d=1/0;             % d is reset for every approach
    Rot_init(:,:,1,k_count)=eye(3);  % Initial rotation matrix created
    Rot_init(:,:,2,k_count)=eye(3);  % Initial rotation matrix created
    Rot_init(:,:,3,k_count)=eye(3);  % Initial rotation matrix created
    P=P_mp;            % Measured points are reset for every approach
    
    
    if k_count==1    % Initial alignment error is not changed if this is the first approach
        T_rand(:,k_count)=[0;0;0];
        
        P=arg.InitialTrans*[P;ones(1,N)];
        P=P(1:3,:);
    else
        
        if k_count>=50 && k_count<=(50+47)
            P=FuncICP_rotMat(pi,'x','r')*P;
            Rot_init(:,:,k_count)=Rot_init(:,:,k_count)*FuncICP_rotMat(pi,'x','r');
        elseif k_count>=98 && k_count<=(98+47)
            P=FuncICP_rotMat(pi,'y','r')*P;
            Rot_init(:,:,k_count)=Rot_init(:,:,k_count)*FuncICP_rotMat(pi,'x','r');
        elseif k_count>(98+47)
            P=FuncICP_rotMat(pi,'z','r')*P;
            Rot_init(:,:,k_count)=Rot_init(:,:,k_count)*FuncICP_rotMat(pi,'x','r');
        end
        
        curRot = rot_map(k_count-1,:);
        
        theta_x=curRot(1);
        R_x=FuncICP_rotMat(theta_x,'x','d');
        
        theta_y=curRot(2);
        R_y=FuncICP_rotMat(theta_y,'y','d');
        
        theta_z=curRot(3);
        R_z=FuncICP_rotMat(theta_z,'z','d');
        
        T_rand(:,k_count)=[(rand-0.5)*L_scale;(rand-0.5)*L_scale;(rand-0.5)*L_scale];
        %T_rand(:,k_count)=[(0.5)*L_scale;(0.5)*L_scale;-(0.5)*L_scale];
        
        P=P-[mu_p_mp(1)*ones(1,N);mu_p_mp(2)*ones(1,N);mu_p_mp(3)*ones(1,N)];
        
        %         if k_to_t==1
        %             P=R_x*P;
        %             Rot_init(:,:,k_count)=Rot_init(:,:,k_count)*R_x;
        %             k_to_t=k_to_t+1;
        %         elseif k_to_t==2
        %             P=R_y*P;
        %             Rot_init(:,:,k_count)=Rot_init(:,:,k_count)*R_y;
        %             k_to_t=k_to_t+1;
        %         else
        %             P=R_z*P;
        %             Rot_init(:,:,k_count)=Rot_init(:,:,k_count)*R_z;
        %             k_to_t=1;
        %         end
        Rot_init(:,:,1,k_count)=Rot_init(:,:,1,k_count)*R_x;
        Rot_init(:,:,2,k_count)=Rot_init(:,:,2,k_count)*R_y;
        Rot_init(:,:,3,k_count)=Rot_init(:,:,3,k_count)*R_z;
        P=R_x*P;
        P=R_y*P;
        P=R_z*P;
        
%         res=22;
%         x = P(1,:);y = P(2,:);z = P(3,:);
%         x = reshape(x',[res,res]);y = reshape(y',[res,res]);z = reshape(z',[res,res]);
%         f2=figure;
%         set(f2,'Color',[1,1,1])
%         hold on;
%         plot3(PP(1,:),PP(2,:),PP(3,:),'.','MarkerSize',0.1)
%         surf(x,y,z);
        P=P+repmat([T_rand(1,k_count)+mu_p_surf(1);T_rand(2,k_count)+mu_p_surf(2);T_rand(3,k_count)+mu_p_surf(3)],1,N);
    end
    
    while converge<2 && d(end)>arg.MinE && k<arg.k_max
        tic;
        k=k+1;
        
        
        % Calculate the corresponding closest points
        [ X,~,~] = FuncICP_closestPoint_m(P_surf,P);
        
        
        % Calculate centres of gravity for both sets of points
        mu_p=sum(P,2)/size(P,2);
        mu_x=sum(X,2)/size(X,2);
        
        % Calculate the cross-covariance matrix
        P_ = P - mu_p;
        X_ = X - mu_x;
        sigma = P_*transpose(X_);
        
        switch arg.Minimize
            case 'SVD'
                [U,~,V]=svd(sigma); % Calculate SVD of the cross-covariance matrix
                qr_ = V*diag([1 1 det(U*V')])*U'; % Compute the optimal rotation matrix
            case 'Quaternion'
                [qr_] = Quat(sigma);
        end
        qr(:,:,k)=qr_;
        
        % Compute the optimal translation vector
        qt(:,k)=mu_x-qr(:,:,k)*mu_p;
        
        % Transform the measured point set
        P=qr(:,:,k)*P + repmat(qt(:,k),1,size(P,2));
        
        % Calculate the mean error
        d(k)=sum(sum((P-X).^2).^0.5)/N;
        
        % Convergence check
        if k==1
            e=d(k);
        else
            e=abs(d(k)-d(k-1))  ;
        end
        k_conv=k;
        if e<arg.ConvergeCrit
            converge=converge+1 ;
            if k_conv==(k_next+1) && converge==2
                terminateCondi_ = 2; % Iteration Termination condition: convergence reached.
            elseif converge==2
                converge=0;
                k_next=k_conv;
            else
                k_next=k_conv;
            end
            
        end
        
        timeCost(k) = toc;
        if k/arg.k_max>nextProgressReport;
            disp(['k count: ',num2str(k),'/',num2str(arg.k_max),...
                '; Approach: ',num2str(k_count),'/',num2str(arg.NoApproach),...
                '; deviation = ',num2str(round(d(k),6))]);
            timeEst = mean(timeCost)*(arg.k_max-k+(arg.NoApproach-k_count)*arg.k_max)/60;
            disp(['Estimated time remaining: ',num2str(round(timeEst,2)),'min']);
            nextProgressReport = nextProgressReport+0.33;
        end
    end
    
    % Iteration Termination condition: minimum error reached
    if d(end)<= arg.MinE
        terminateCondi_ = 1; 
        disp(['Minimum error reached!',num2str(d(end))]);
    elseif terminateCondi_~= 2 % Iteration Termination condition: if not converged
        terminateCondi_ = 3; % Iteration Termination condition: maximum iteration reached
    else
        disp(['Convergence reached!',num2str(d(end))]);
    end
    
    
    k_col(k_count)=k;       % Number of iterations of the approach collected
    % Transolation matrix of the approach collected
    Trans_collect{k_count}=qt(:,:);
    % Rotational array of the approach collected
    Rotation_collect{k_count}=qr(:,:,:);
    % The final mean error of the apporach collected
    d_col(k_count)=d(end);
    icp.d{k_count}=d;
    icp.terminateCondi{k_count} = terminationInfo{terminateCondi_};
    icp.timeCost_col{k_count}=timeCost;
    % Parameters for the next approach set
    k_count=k_count+1;
    k=0;
    e=1/0;
    converge=0;
    
    
    %     disp(['k count: ',num2str(k_count), '; deviation = ',num2str(round(d_col(k_count),6))]);
    %     timeCost(k_count) = toc;
    %     disp(['Average time cost: ',num2str(round(mean(timeCost),6)),'s']);
end

%% The following section calcualtes the best transformation matrix found
% The best aproach is found by selecting the set corresponding to the
% minimum mean error
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
    trans(:,:)=trans(:,:)*arg.InitialTrans;
else
    trans(:,:)=trans(:,:)*[eye(3) (T_rand(:,I)+mu_p_surf);0 0 0 1];
    trans(:,:)=trans(:,:)*[Rot_init(:,:,1,I) zeros(3,1);0 0 0 1];
    trans(:,:)=trans(:,:)*[Rot_init(:,:,2,I) zeros(3,1);0 0 0 1];
    trans(:,:)=trans(:,:)*[Rot_init(:,:,3,I) zeros(3,1);0 0 0 1];
    trans(:,:)=trans(:,:)*[eye(3) mu_p_mp*(-1);0 0 0 1];
end
icp.Trans=Trans_collect;
icp.Rot=Rotation_collect;
end


%%
function [R] = Quat(Sigma)
% The anti-symmetric part of the cross covariance matrix is found
Delta=Sigma-Sigma';
% The cyclic components are found from the anti-symmetric matrix
Delta=[Delta(2,3) Delta(3,1) Delta(1,2)]';
% The Q matrix is constructed
B=Sigma+Sigma'-trace(Sigma)*eye(3);
Q=[trace(Sigma) Delta';Delta B];

% Eigen vectors and eigenvalues are found
[eig_Vec,eig_Val]=eig(Q);

% The maximum eigenvalue is found
C_1=max(eig_Val);
[~,C_2]=max(C_1);
% The unit quaternion vector is defined as the eigenvector
% correspoinding to the maximum eigenvalue
q=eig_Vec(:,C_2);

% The optimum rotational matrix is found
R=[q(1)^2+q(2)^2-q(3)^2-q(4)^2 2*(q(2)*q(3)-q(1)*q(4)) 2*(q(2)*q(4)+q(1)*q(3));
    2*(q(2)*q(3)+q(1)*q(4)) q(1)^2-q(2)^2+q(3)^2-q(4)^2 2*(q(3)*q(4)-q(1)*q(2));
    2*(q(2)*q(4)-q(1)*q(3)) 2*(q(3)*q(4)+q(1)*q(2)) q(1)^2-q(2)^2-q(3)^2+q(4)^2];
end