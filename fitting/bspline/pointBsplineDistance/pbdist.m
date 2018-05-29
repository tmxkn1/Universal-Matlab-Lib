function [pdist, pclst] = pbdist(p, c, knotu, knotw, ku, kw, varargin)
%PBDIST Evaluation of distances between points and a B-spline surface
% PDIST = PBDIST(P, C, KNOTU, KNOTW, KU, KW) returns distances PDIST 
% between each point in P (x-by-3) and a B-spline surface defined by 
% control points C (n-by-m-3), with knots KNOTU (n+ku-by-1) and of order 
% KU in u-direction and knots KNOTW (m+kw-by-1) and of order KW in 
% w-direction.
%
% [PDIST, PCLST] = PBDIST(...) also returns the closest points on the
% B-spline surface.
%
% The function first decomposes the B-spline surface into several Bezier
% patches and generate a number of points on each patch. These generated
% points are compared with P until the closest point for each point in P 
% found. This establishes initial guesses of projections of P on the
% B-spline surface.
% Next, the initial guesses are refined by subdividing the Bezier patches
% around the closest points until the distances between the closest points
% and P converge.
%
% The number of points generated on EACH Bezier patch can be specified with
% 'StartPtsSize':
% [...] = PBDIST(..., 'StartPtsSize', [A, B])
% A defines the number of points in u-direction and B, w-direction. By
% default, 15-by-30 points are generated.
%
% The converging criteria can be specified with 'ConvCriteria':
% [...] = PBDIST(..., 'ConvCriteria', A)
% By default, A=0.001.
%
% [...] = PBDIST(..., 'ProgressBar', 'on') shows a progress bar.
%
% -------
% Required custom functions:
% bsp2bez, bezsubdiv
%
% -------
% Required Toolboxes:
% Statistics and Machine Learning Toolbox
% Parallel Computing Toolbox (can be removed by changing 'parfor' to 'for'.
%
% ref: 
% Ma, Y. & Hewitt, W. (2003) Point inversion and projection for NURBS 
% curve: Control polygon approach, Theory and Practice of Computer Graphics
% TPCG 2003, vol. 20, pp. 113¨C120.
%
% Coded by Zhengyi Jiang from The University of Manchester

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