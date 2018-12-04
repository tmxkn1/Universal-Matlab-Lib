function [pdist, pclst] = pbdist(p, c, knotu, knotw, ku, kw, varargin)
% Calculates the distances between points in 3D and a B-spline surface.
%
% The function first decomposes the B-spline surface into several Bezier
% patches and generates a number of points (15*30 by default) on each
% patch. Then, for each point in the input point set, a closest point is
% found on the Bezier patch. Next, for each closest point, the Bezier patch
% it belongs to is subdivided around the point and a new closest point is
% found for the corresponding point in the point set. This is performed
% iteratively until the distance the two points converge.
%
% -------------------------------------------------------------------------
% USE:
%
% D = pbdist(P, C, VU, VW, KU, KW) returns distances, D, between each point 
% in P and a B-spline surface defined by control points, C, with knot
% vectors, KU and KW, in u- and w-direction.
%
%   Input:  P  - a p*3 matrix, containing a point set.
%           C  - an n*m matrix, containing the control points.
%           VU - a column vector, containing the parametrised coordinates
%                on the B-spline surface in the u-direction.
%           VW - a column vector, containing the parametrised coordinates
%                on the B-spline surface in the w-direction.
%           KU - a column vector, containing the knots of the B-spline 
%                surface in the u-direction.
%           KW - a column vector, containing the knots of the B-spline 
%                surface in the w-direction.
%
% [D, PC] = pbdist(...) also returns the coordinates of the corresponding 
% closest points, PC, on the B-spline for surface for each point in the 
% point set.
%
% The number of points generated on EACH Bezier patch can be specified with
% 'StartPtsSize':
% [...] = pbdist(..., 'StartPtsSize', [A, B])
% A represents the number of points in the u-direction and B represents the
% number in the w-direction.
%
% The converging criteria can be specified with 'ConvCriteria':
% [...] = pbdist(..., 'ConvCriteria', C)
% By default, C = 0.001.
%
% [...] = pbdist(..., 'ProgressBar', 'on') also shows a progress bar.
%
% -------------------------------------------------------------------------
% Required custom functions:
% bsp2bez, bezsubdiv
%
% Toolboxes:
% Statistics and Machine Learning Toolbox
% Parallel Computing Toolbox (can be removed by changing 'parfor' to 'for'.
%
% -------------------------------------------------------------------------
% Reference: 
%
% Ma, Y. & Hewitt, W. (2003) Point inversion and projection for NURBS 
% curve: Control polygon approach, Theory and Practice of Computer Graphics
% TPCG 2003, vol. 20, pp. 113:C120.

pars = inputParser;
addRequired(pars,'p');
addRequired(pars,'c');
addRequired(pars,'knotu');
addRequired(pars,'knotw');
addRequired(pars,'ku');
addRequired(pars,'kw');
addParameter(pars,'StartPtsSize',[15,30]);
%addParameter(pars,'SubdividePtsSize',[15,30],@(x) numel(x)==2); % unused.
addParameter(pars,'ConvCriteria', 0.001);
parse(pars, p, c, knotu, knotw, ku, kw, varargin{:});
arg = pars.Results;

convcrit = arg.ConvCriteria;
pchsize = arg.StartPtsSize;
%spchsize = arg.SubdividePtsSize;

%% Check everything before we start
if ~isnumeric(p)
    error('Point set P must be a n-by-3 or 3-by-n matrix.');
elseif size(p,2) ~= 3
    if size(p,1) == 3
        p = p';
    else
        error('Point set P must be a n-by-3 or 3-by-n matrix.');
    end
end
if ~(size(c,3) == 3 && isnumeric(c))
    error('Control points C must be a n-by-m-by-3 3-D matrix.');
end
if ~(numel(knotu) > 1 && isnumeric(knotu))
    error('Knot vector KNOTU must be a vector.');
end
if ~(numel(knotw) > 1 && isnumeric(knotw))
    error('Knot vector KNOTW must be a vector.');
end
if ~(isscalar(ku) && isnumeric(ku))
    error('Order of the surface KU must be a numerical scalar.');
end
if ~(isscalar(kw) && isnumeric(kw))
    error('Order of the surface KW must be a numerical scalar.');
end
if ~(numel(pchsize)==2 && isnumeric(pchsize))
    error('Argument for ''StartPtsSize'' must be a 1-by-2 numerical vector.');
end
if ~(isscalar(convcrit) && isnumeric(convcrit))
    error('Argument for ''ConvCriteria'' must be a scalar.');
end
 
%% Now begin
% number of points in a patch 
%ptnum = pchsize(1)*pchsize(2);

% u/w offset
if convcrit > 0.001
    uwoffset = 0.0001;
else
    uwoffset = convcrit * 0.1;
end

% b-spline to bezier
[bezctrlpt] = bsp2bez(knotu, knotw, c, ku, kw);
bezctrlpt = permute(bezctrlpt,[2,1,3]);
pchnumu = size(bezctrlpt,1);
pchnumu = (pchnumu-1)/(ku-1);
%pchnumw = (pchnumw-1)/(kw-1);

% fit Bezier patches
uvals = linspace(0+uwoffset/2,1-uwoffset/2,pchsize(1));
wvals = linspace(0+uwoffset/2,1-uwoffset/2,pchsize(2));
[bezx,bezy,bezz] = bezfit(bezctrlpt,ku,kw,...
    uvals, wvals);
bezp = [bezx(:),bezy(:),bezz(:)];
% search for nearest points
kdOBJ = KDTreeSearcher(bezp);
[matchlist, distlist] = knnsearch(kdOBJ,p);
forcequitc = 20;

% pre-allocation
pdist = zeros(size(p,1),1);
pclst = zeros(size(p));

% total number of points
np = numel(matchlist);

%% go through each point
parfor i = 1:np
    % get the next point for the loop
    m = matchlist(i); d = distlist(i); p_ = p(i,:);
    % reset some variables
    cpmat = bezctrlpt; bp = bezp;
    pchs = pchsize; pchnu = pchnumu; %spchs = spchsize;
    uvs = uvals; wvs = wvals; 
    %{     
    figure; hold on;axis equal
    plot3(p_(:,1),p_(:,2),p_(:,3),'o','Color',[1,0,0],'MarkerFaceColor',[0,0,1]);
    plot3(bezp(m,1), bezp(m,2),bezp(m,3),'^k','Color',[1,0,0],'MarkerFaceColor',[0,0,1]);
    surf(bezx,bezy,bezz)
    %}
    % preset some variables
    d_= d; mp_= bp(m,:);
    % initialise the loop
    ic=1; err = 1;
    while err>convcrit       
        % find out the location of the point in the entire Bezier surface
        pcolid = ceil(m/pchs(1)/pchnu); % column id of the point
        prowid = m-(pcolid-1)*pchs(1)*pchnu; % row id of the point
        % and the patch
        pchrowid = ceil(prowid/pchs(1)); % row id of the patch
        pchcolid = ceil(pcolid/pchs(2)); % column id of the patch
        % find the location of the point in the patch
        uc = prowid-(pchrowid-1)*pchs(1);
        wc = pcolid-(pchcolid-1)*pchs(2);
        % extract all neighbouring points
        uc = uc-1:uc+1; wc = wc-1:wc+1;
        if uc(1) == 0; uc(1) = []; end
        if uc(end) > numel(uvs); uc(end) = []; end
        if wc(1) == 0; wc(1) = []; end        
        if wc(end) > numel(wvs); wc(end) = []; end
        u = uvs(uc); w = wvs(wc);
        % extract the control points if 
        if ic == 1
            cpmat = cpmat((pchrowid-1)*(ku-1)+1:pchrowid*(ku-1)+1,...
                (pchcolid-1)*(kw-1)+1:pchcolid*(kw-1)+1,:);
        end
        % subdivide at minimum u and minimum w, take the bottom right patch
        [~,~,~,~,cp] = bezsubdiv(cpmat,min(u),min(w));
        % subdivide the patch at max u and max w, take the top left patch
        [~,cpmat] = bezsubdiv(cp,(max(u)-min(u))/(1-min(u)),(max(w)-min(w))/(1-min(w)));
        % grow point number per patch
        pchs = ceil(pchs*1.1);
        uvs = linspace(0+uwoffset/2,1-uwoffset/2,pchs(1));
        wvs = linspace(0+uwoffset/2,1-uwoffset/2,pchs(2));
        [bx,by,bz] = bezfit(cpmat,ku,kw,uvs, wvs);
        bp = [bx(:),by(:),bz(:)];
        % search for nearest points
        kdOBJ = KDTreeSearcher(bp);
        [m, d] = knnsearch(kdOBJ,p_);
        %{
        surf(cp(:,:,1),cp(:,:,2),cp(:,:,3))
        surf(cpmat(:,:,1),cpmat(:,:,2),cpmat(:,:,3))
        surf(bx,by,bz,'EdgeColor','[1,1,1]','FaceColor',[0.7,0.7,0.7]);
        plot3(bp(m,1),bp(m,2),bp(m,3),'^k','Color',[1,0,0],'MarkerFaceColor',[0,0,1])
        %}
        err = abs(d_(ic)-d);
        ic = ic+1;
        d_(ic) = d; mp_(ic,:) = bp(m,:);
        if ic>forcequitc
            break;
        end
        pchnu = 1;
    end
    % get the minimum distance, hence the closest point
    pdist(i) = min(d_);
    mp_ = mp_(pdist(i)==d_,:);
    pclst(i,:) = mp_(1,:);
end