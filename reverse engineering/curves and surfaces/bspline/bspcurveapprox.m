function [x,y,z,knotV, ctrlp] = bspcurveapprox (x,y,z, varargin)
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

p = bspfit(ufit, ctrlp, k, knotV);
x = p(:,1);
y = p(:,2);
z = p(:,3);

