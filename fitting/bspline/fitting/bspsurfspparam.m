function [uu, ww, p_ordered, tri] = bspsurfspparam(u,w,p,plot2dmesh)
% Calculates Floater's shape preserving parametrisation.
%
% Floater's shape preserving parametrisation requires a surface mesh of
% the 3D points. This function does not require you to supply the surface
% mesh but instead a set of coarsely parametrised coordinates using other
% methods such as chord-length parametrisation.
%
% -------------------------------------------------------------------------
% USE:
%
% [UU, WW] = bspsurfspparam(U, W, P) parametrises a set of 3D points P 
% using Floater's shape preserving parametrisation algorithm. U contains 
% the coarsely parametrised coordinates in the u-direction and W contains 
% those in the w-direction. The output is a set of shape-preserved 
% parametrised coordinates, UU and WW, in the u- and w-direction
% respectively.
%
% -------------------------------------------------------------------------
% See also: bspsurffit, bspcurvefit

if nargin < 4
    plot2dmesh = 0;
end

%% Triangulation
% combine u and w
uw = [u, w];
uw_inside = uw;
p_inside = p;

% find boundary
bu1i = find(uw_inside(:,2)==0);
bu1_ = uw_inside(bu1i,:);
uw_inside(bu1i,:) = [];
pu1_ = p_inside(bu1i,:);
p_inside(bu1i,:) = [];

bu2i = find(uw_inside(:,2)==1);
bu2_ = uw_inside(bu2i,:);
uw_inside(bu2i,:) = [];
pu2_ = p_inside(bu2i,:);
p_inside(bu2i,:) = [];

bw1i = find(uw_inside(:,1)==0);
bw1_ = uw_inside(bw1i,:);
uw_inside(bw1i,:) = [];
pw1_ = p_inside(bw1i,:);
p_inside(bw1i,:) = [];

bw2i = find(uw_inside(:,1)==1);
bw2_ = uw_inside(bw2i,:);
uw_inside(bw2i,:) = [];
pw2_ = p_inside(bw2i,:);
p_inside(bw2i,:) = [];

u_boundary=[bu1_;bw1_;bu2_;bw2_];
uw_c = [uw_inside;u_boundary];
tri=delaunay(uw_c);

if plot2dmesh
    figure('Color',[1,1,1]);hold on;
    set(gcf,'Position',[0,460,315,309]);
    trimesh(tri,uw_c(:,1),uw_c(:,2),'Color',[0.3,0.3,0.3]);              
    %highlight first line
    I = find(uw_c(:,1) == 0);
    plot(uw_c(I,1),uw_c(I,2),'+-','Color',[162/255,20/255,47/255],...
        'LineWidth',2)
    
    axis equal;axis tight;box on;grid on;
    xlabel x;ylabel y;
    set(gca,'LooseInset',get(gca,'TightInset'));
    
end

p_ordered = [p_inside;pu1_;pw1_;pu2_;pw2_];

%% Shape-preserving parameterization
n=size(u_boundary,1);
x = p_ordered(:,1);
y = p_ordered(:,2);
z = p_ordered(:,3);

[tp,td]=bspsurfspparamedge(tri,x,y,z,n);   % compute barycentric coordinates of data points

pu=tp\(td*u_boundary(:,1));      % compute the u and w value for each point
pw=tp\(td*u_boundary(:,2));

puw=[pu,pw; u_boundary];

uu=puw(:,1);
ww=puw(:,2);

%%
function[tp,td]=bspsurfspparamedge(t,x,y,z,np1)  
% compute the barycentric coordinate of data points

%% Sorting the edges of each polygon
np=size(x,1);
N=np-np1;
v=ones(1,N);
tp=diag(v);
td=zeros(N,np1);   %% Initialize the variables

n=1;
for id=1:N   % Check every point in sequence
    edge=[]; ve=[]; angle=[];
    lin=[]; newve=[]; x1=[]; y1=[]; z1=[];
    
    %% STEP 1. Define the polygon
    ppx=0;
    ppy=0;
    angleplane=0;
    % Initialize variables at start of each loop
    checkt=any(t==id,2);  %% extracting the polygon of which the center is "id" tri=t(checkt,:);
    tri=t(checkt,:);
    nt=size(tri,1);
    
    if isempty(tri)%% In case of noisy point existing, the program will skip the noisy point
        continue
    end
    % figure(8)
    % trimesh(tri,x,y,z)
    checkt=true(nt,1);
    
    for j=1:nt        %% position the center point of the polygon in the first column
        if tri(j,2)==id;
            tri(j,2)=tri(j,1);
        elseif tri(j,3)==id;
            tri(j,3)=tri(j,1);
        end
        tri(j,1)=id;
    end
    
    % extracting all edges of one polygon
    e=[tri(:,[1,2]); tri(:,[1,3])];
    l = 1:nt;
    [e,j,j] = unique(e,'rows');  %% Remove duplicated edges
    e2 = [j(l), j(l+nt)];   % e2 stores the relationship between the triangle and its edges
    ne=size(e,1);
    e3 = zeros(ne,2);
    count= ones(ne,1);
    
    % e3 defines the two triangles owned by one edge
    for l=1:nt
        for j=1:2
            ce = e2(l,j);
            e3(ce,count(ce))=l;
            count(ce)=count(ce)+1;
        end
    end
    
    % position center point in the first column
    for j=1:ne     
        if e(j,2)==id
            e(j,2)=e(j,1);
        end
        e(j,1)=id;
    end
    k=1;
    m=1;                    % m is the ID of the triangle that will be extracted
    n=1;
    edge(k,:)=e(e2(m,1),:); % define the first edge that is selected
    idp(k)=edge(k,2);       % take all boundary points of the polygon
    k=k+1;
    
    while any(checkt==1)
        tri1=e3(e2(m,n),1);    % Select two adjacent triangles which share a common edge    tri2=e3(e2(m,n),2);
        tri2=e3(e2(m,n),2);
        checkt(m)=false;
        if  tri1==m             % If tri1 is the last triangle
            ed1=e(e2(tri2,1),:);
            ed2=e(e2(tri2,2),:); % one edge of tri2 is the next edge in sequence
            if ed1(2)==edge(k-1,2)
                edge(k,:)=ed2;
                idp(k)=ed2(2);
                k=k+1;
                m=tri2;
                n=2;
            elseif ed2(2)==edge(k-1,2)
                edge(k,:)=ed1;
                idp(k)=ed1(2);
                k=k+1;
                m=tri2;
                n=1;
            end
        else              % If tri2 is the last triangle
            ed1=e(e2(tri1,1),:);
            ed2=e(e2(tri1,2),:);  % one edge of tri2 is the next edge in sequence
            if ed1(2)==edge(k-1,2)
                edge(k,:)=ed2;
                idp(k)=ed2(2);
                k=k+1;
                m=tri1;
                n=2;
                
            elseif ed2(2)==edge(k-1,2)
                edge(k,:)=ed1;
                idp(k)=ed1(2);
                k=k+1;
                m=tri1;
                n=1;
                
            end
        end
        
    end
    k=k-1;
    
    %% STEP 2. Mapping the space polygon into plane
    % Extract the coordinate of the end points of edges
    for j=1:k
        x1(j,:)=x(edge(j,:));  
        y1(j,:)=y(edge(j,:));
        z1(j,:)=z(edge(j,:));
    end
    
    % Calculate the length of each eage
    for j=1:k                 
        ve(j,:)=[x1(j,2)-x1(j,1),y1(j,2)-y1(j,1),z1(j,2)-z1(j,1)];
        lin(j)=norm(ve(j,:));
    end
    
    % Calculate the sum of the solid angles of two adjacent edge
    sumangle=0;
    for j=1:k-1              
        result = dot(ve(j+1,:),ve(j,:))/(norm(ve(j+1,:))*norm(ve(j,:)));
        if result>1 % control the accuracy of dot
            result = 1;
        end
        angle(j)=acos(result);
        sumangle=angle(j)+sumangle;  % plus the solid angle
    end
    
    % Solve Eq.1 of shape-preserving parameterization
    for j=1:k-1  
        angleplane(j+1)=angle(j)/sumangle*2*pi;
    end
    newve(1,:)=[1,0]; % reference vector
    nor=1;
    
    % Calculate transformation matrix
    for j=1:k                 
        cs=cos(angleplane(j));
        si=sin(angleplane(j));
        
        ro=[cs,si;-si, cs];
        
        newve(j+1,:)= newve(j,:)/nor*ro.*lin(j); % map the polygon into plane
        nor=lin(j);  
    end
    
    % Store the coordinate of vertices of polygon
    for j=2:k
        ppx(j,:)=newve(j,1)+ppx(1);    
        ppy(j,:)=newve(j,2)+ppy(1);
        pp=[ppx,ppy];
    end
    %figure(5)
    %plot(ppx,ppy,'*-');
    
    %% STEP 3. Get final stuff done
    % Calculate the Barycentric Coordinates
    [tarea]=bspsurfspparamarea(ppx,ppy); 
    for l=1:k-1                  % record �� into tp and td
        if idp(l)<=N
            tp(id,idp(l))=tp(id,idp(l))-tarea(l+1);
        else
            idp(l)=idp(l)-N;
            td(id,idp(l))=tarea(l+1);
        end
        
    end
end

function [tarea]=bspsurfspparamarea(ppx,ppy)

ne=size(ppx,1);   
n11=0;
n22=0;
t=zeros(ne,1);
ppx(ne+1:2*ne-1)=ppx(2:ne);   %% order the vertices of polygon
ppy(ne+1:2*ne-1)=ppy(2:ne);
for i=2:2*ne-1
ve(i,:)=[ppx(i)-ppx(1),ppy(i)-ppy(1)]; %% calculate the length from vertices to center point
length(i)=norm(ve(i,:));
end

%Let the first vertice selected be the 2nd one in the list of the polygon
%points. then the first i=4 and i-1=3, first k=2, therefore, the first set
%of points are 2,3,4.
o=4;
l=ne;

for k=2:ne
for i=o:l
%Select three consecutive vertices around the centre point, x, and calculate the
%area of the triangle formed by these three vertices. 
  ve1=[ppx(i-1)-ppx(k),ppy(i-1)-ppy(k)]; 
  length1=(norm(ve1));
  
  ve2=[ppx(i)-ppx(k),ppy(i)-ppy(k)];  
  length2=(norm(ve2)) ;
  
  ve3=[ppx(i)-ppx(i-1),ppy(i)-ppy(i-1)];
  length3=(norm(ve3));
  
  %Calculate the area of triangle formed by two vertices and the centre
  %point, devided by the area of the triangle calculated in the previous
  %step.
  p=(length1+length2+length3)/2;
  s=sqrt(p*(p-length1)*(p-length2)*(p-length3));% calculate the area of the triangle 
  
  p=(length(i)+length(i-1)+length3)/2;
  s1=sqrt(p*(p-length3)*(p-length(i))*(p-length(i-1)));
  t1=abs(s1/s);                              % t1: the barycentric coordinates of x
  
  p=(length2+length(i)+length(k))/2;
  s2=sqrt(p*(p-length(k))*(p-length(i))*(p-length2));
  t2=abs(s2/s);                              % t2: the barycentric coordinates of x
  
  p=(length1+length(i-1)+length(k))/2;
  s3=sqrt(p*(p-length1)*(p-length(i-1))*(p-length(k)));
  t3=abs(s3/s);                              % t3: the barycentric coordinates of x
  
  %Add up those results, i.e. t1, t2, t3, if the sum equals 1, it means the
  %centre point is inside the triangle formed by the three vertices
  %selected in the first step. And t1, t2 and t3 are therefore barycentric
  %coordinates of the centre point.
  tsum=t1+t2+t3;      
  if abs(tsum-1)<0.00001   %% if the x located in triangle,"tsum"will equal to 1 
      if i==ne+1
          n1=2;
          n2=ne;
      elseif i>ne+1
          n1=i-ne+1;
          n2=i-ne;
      else
          n1=i;
          n2=i-1;
      end
     
      t(k)=t(k)+t1;
      t(n2)=t(n2)+t2;
      t(n1)=t(n1)+t3;
     
    break     
  end 
end
o=o+1;
l=l+1;
end
tarea=t/(ne-1);