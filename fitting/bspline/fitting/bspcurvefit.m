function [x,y,z,knotV, ctrlp] = bspcurvefit (x,y,z, varargin)
% Fits a set of 3D points with a B-spline curve.
%
% -------------------------------------------------------------------------
% USE:
%
% [XO, YO, ZO] = bspcurvefit(X, Y, Z) fit a B-spline curve to a set of 
% points whose coordinates are X, Y and Z. The coordinates of the fitted 
% B-spline curve are XO, YO and ZO. 
%
% The number of points on the B-spline curve can be defined by supplying a 
% parameter 'NumOfNewPoints' with a value:
%       [...] = bspcurvefit(..., 'NumOfNewPoints', n);
% If unspecified, the number of points on the B-spline curve is equal to
% twice the number in the original point set.
%
% The order of the B-spline curve can be defined by supplying a parameter
% 'Order' with a value:
%       [...] = bspcurvefit(..., 'Order', k);
% If unspecified, the order of the curve is four.
%
% The number of segments on the B-spline curve can be defined by supplying
% a parameter 'NumOfSegments' with a value:
%       [...] = bspcurvefit(..., 'NumOfSegments', s);
%
% [..., V] = bspcurvefit(...) also outputs the knot vector of the B-spline
% curve.
%
% [..., C] = bspcurvefit(...) also outputs the control points used in the 
% B-spline fitting.
%
% -------------------------------------------------------------------------
% See also: bspbasis, bspchordlparam, bspcontrolpt, bspgenpoints,
% bspknot, bspsurffit, bezfit

pars = inputParser;
addRequired(pars,'x');
addRequired(pars,'y');
addRequired(pars,'z');
addParameter(pars,'NumOfNewPoints',0);
addParameter(pars,'Order',4);
addParameter(pars,'NumOfSegments',0);
parse(pars,x,y,z,varargin{:});
arg = pars.Results;

nou = arg.NumOfNewPoints;
k = arg.Order;
noseg = arg.NumOfSegments;
if nou == 0
    nou = numel(x)*2;
end
if noseg == 0
    noseg = ceil(numel(x)/3) - k;
end

noctrlp = noseg + k - 1;
u = bspchordlparam([x,y,z], noseg);
[knotV] = bspknot(noctrlp, k);
ctrlp = bspcontrolpt([x,y,z], u, noctrlp, k, knotV);

ufit = linspace(0, noseg, nou);

p = bspgenpoints(ufit, ctrlp, k, knotV);
x = p(:,1);
y = p(:,2);
z = p(:,3);