function fh = plotpbdist(fh, dist, p, pmatched)
% Plots the distances from one set of points to another set of points
%
% This function does not calculate the distances. 
% -------------------------------------------------------------------------
% USE:
%
% H = plotpbdist(D, P, PM) plots a colourful figure shows the distances, D,
% from a set of points, P, to another set of points, PM, where 
%       PM(x,y,z)_i = P(x,y,z)_i + D*n(x,y,z)_i.
% The output is a handle to the figure object.
%
% H = plotpbdist(H, D, P, PM) plots the figure on an existing figure object
% whose handle is H.

if nargin ~= 4
    pmatched = p;
    p = dist;
    dist = fh;
    fh = figure;
end
hold on;
cm = jet(101);
mind = min(dist);
maxd = max(dist);
d = maxd-mind;

xc = p(:,1);yc = p(:,2);zc = p(:,3);
xp = pmatched(:,1);yp = pmatched(:,2);zp = pmatched(:,3);
x_ = [xc';xp';nan(size(xc'))];
y_ = [yc';yp';nan(size(xc'))];
z_ = [zc';zp';nan(size(xc'))];

for i = 1:numel(dist)
    c = cm(round((dist(i)-mind)/d*100)+1,:);
    plot3(x_(:,i),y_(:,i),z_(:,i),'.-','Color',c)
end
colormap(cm);
colorbar;
% colorbar limits 
caxis([mind, maxd])
% ylabel(hc,...
%     sprintf('(mm) RMSD: %.3f; Max: %.3f; Min: %.3f; Mean: %.3f; STD: %.3f',...
%     sqrt(mean((dist).^2)),maxd,mind,mean(dist),std(dist)))
axis equal; axis tight; hold off;