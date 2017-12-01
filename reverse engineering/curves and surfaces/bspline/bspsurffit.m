function [varargout]=bspsurffit(varargin)
% BSPSURFFIT calculates a B-spline surface based on a set of input
% parameters.
%
% [BS] = BSPSURFFIT(P,U,W,CU,CW) 
% requires inputs:
%   * P    a set of surface points in 3D Cartesian coordinates
%   * U    parameters in u-direction
%   * W    parameters in w-direction
%   * CU   Number of control points in u-direction
%   * CW   Number of control points in w-direction
% returns a structured varible BS:
%   * BS.cx, BS.cy, BS.cz are CU-by-CW matrices for the x, y and z 
%     components of the control points, respectively.
%   * BS.knotu and BS.knotw are knot vectors in u and w directions,
%   respectively.
%
% [FP, BS] = BSPSURFFIT(...,NU,NW) 
% with two more parameters:
%   * NU   Number of fitted points in u direction
%   * NW   Number of fitted points in w direction
% returns fitted surface points, FP (NU*NW-by-3), as the first argument,
% and the same structured variable metioned above.
%
% by Zhengyi Jiang, The University of Manchester, 2017

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