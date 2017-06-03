% In this demo, a bezier patch of order 4 in u-direction and order 3 in
% w-direction is divided into four bezier patches of the same orders.

load('bezier patch')

% split the Bezier patch at u=0.5 and w=0.5
u = 0.5; w = 0.5;
[CP,CLU,CRU,CLW,CRW] = bezsubdiv(c,u,w);

% plot
figure; hold on; 
mesh(c(:,:,1),c(:,:,2),c(:,:,3),'Marker','+','FaceColor','none',...
    'EdgeColor',[0,0.45,0.74]);
surf(CLU(:,:,1),CLU(:,:,2),CLU(:,:,3),'FaceColor',[1 0.84 0])
surf(CRU(:,:,1),CRU(:,:,2),CRU(:,:,3),'FaceColor',[0.93 0.69 0.13])
surf(CLW(:,:,1),CLW(:,:,2),CLW(:,:,3),'FaceColor',[0.87 0.49 0])
surf(CRW(:,:,1),CRW(:,:,2),CRW(:,:,3),'FaceColor',[0.85 0.33 0.1])
mesh(CP(:,:,1),CP(:,:,2),CP(:,:,3),'Marker','o','FaceColor','none',...
    'EdgeColor',[0.2,0.2,0.2]);
legend('Original control polygon','Divided control polygon',...
    'Divided: Upper left','Divided: Upper right','Divided: bottom left',...
    'Divided: bottom right');
xlabel x; ylabel y; zlabel z;
axis equal; grid on; axis tight;
view(-57,64)
    