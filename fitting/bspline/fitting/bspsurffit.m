function [varargout]=bspsurffit(varargin)
% Fits a B-spline surface to a set of 3D points.
%
% -------------------------------------------------------------------------
% USE:
%
% [BS] = bspsurffit(P, U, W, CU, CW) 
%
%   Input:  P  - a set of 3D points.
%           U  - parametric coordinates of P in the u-direction.
%           W  - parametric coordinates of P in the w-direction
%           CU - number of control points in the u-direction.
%           CW - number of control points in the w-direction.
% 
%   Output: BS - a structure with following fields:
%       * cx - CU*CW matrices, containing the X coordinates of the control
%               points.
%       * cy - CU*CW matrices, containing the Y coordinates of the control
%               points.
%       * cz - CU*CW matrices, containing the Z coordinates of the control
%               points.
%       * knotu - knots of the B-spline surface in the u-direction.
%       * knotw - knots of the B-spline surface in the w-direction.
%
% [PG, BS] = bspsurffit(..., NU, NW) by supplying two extra values:
%   NU, the number of points to be generated in the u-direction, and 
%   NW, the number of points to be generated in the w-direction,
% a set of points, PG, can be generated on the B-spline surface and 
% returned.
%
% -------------------------------------------------------------------------
% See also: bspsurfspparam, bspsurfspparam, bspcurvefit


% Parse arguments
if nargin ~= 7 && nargin ~= 9
    error('Not enough input arguments.');
end

if nargout > 2
    error('Too many output arguments.');
end

QX = varargin{1}(:,1);
QY = varargin{1}(:,2);
QZ = varargin{1}(:,3);
u = varargin{2};
w = varargin{3};
ncpu = varargin{4};
ncpw = varargin{5};
ku = varargin{6};
kw = varargin{7};

%% U direction
nosu = ncpu-ku+1; % number of segments

% make sure the last u value equals the number of segments
umax = max(u);
if umax ~= nosu
    u = u / umax * nosu;
end

knotVu = bspknot(ncpu, ku);
U = bspbasis(u, ncpu, ku, knotVu);

%% W direction
nosw = ncpw-kw+1; % number of segments

% make sure the last w value equals the number of segments
wmax = max(w);
if wmax ~= nosw
    w = w / wmax * nosw;
end

[knotVw] = bspknot(ncpw, kw);
W = bspbasis(w, ncpw, kw, knotVw);
W = W';

%% Control points
% repeat each U value by number of control points in w in dim2
UU = kron(U,ones(1,ncpw));
% repeat each W value by number of control poitns in w in dim1
WW = repmat(W,ncpu,1)';

UW = UU.*WW;

% Calculate the control points
CX=UW\QX;
CY=UW\QY;
CZ=UW\QZ;

% Rearrange the control points into matrix form
CXX = reshape(CX,ncpw,ncpu)';
CYY = reshape(CY,ncpw,ncpu)';
CZZ = reshape(CZ,ncpw,ncpu)';

%% Generate new points and output
if nargout == 1 % output structured B-spline surface parameters
    bsp.cx = CXX; bsp.cy = CYY; bsp.cz = CZZ;
    bsp.knotu = knotVu; bsp.knotw = knotVw;
    varargout{1} = bsp;
else % output B-spline surface points matrix and parameters
    % generate uniform u,w for fitting.
    ufit = linspace(0,nosu,varargin{8});
    wfit = linspace(0,nosw,varargin{9});

    Ufit = bspbasis(ufit, ncpu, ku, knotVu);
    Wfit = bspbasis(wfit, ncpw, kw, knotVw);

    FX = Ufit*CXX*Wfit';
    FY = Ufit*CYY*Wfit';
    FZ = Ufit*CZZ*Wfit';

    varargout{1} = [FX(:), FY(:), FZ(:)];
    bsp.x = FX; bsp.y = FY; bsp.z = FZ;
    bsp.cx = CXX; bsp.cy = CYY; bsp.cz = CZZ;
    bsp.knotu = knotVu; bsp.knotw = knotVw;
    varargout{2} = bsp;    
end