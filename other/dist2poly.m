function [d, p] = dist2poly(P,xy)
% DIST2POLY calculates Euclidean distance from a point to a polynomial.
%
% D = dist2poly(P,XY) calculates Euclidean distance, D, from a given point,
% XY, to a polynomical, P(x) defined by the coefficients P.
%
% XY must be a row vector of [X, Y] or a matix of [X0, Y0; X1, Y1 ...]
%
% [D, P] = dist2poly(...) also returns the closest point on the polynomial.

% The distance is calculated by minimizing (x-x0)^2 + (Px-y0)^2 = d^2.
% This can be achieved by differentiating the equation with respect to x:
%        (x-x0) + (Px - y0)*dPx/dx = 0

if size(xy,2) ~= 2
    error('D = dist2poly(P,XY) XY must be a row vector of [X, Y] or a matix of [X0, Y0; X1, Y1 ...]')
end

% fminsearch is faster if degree is greater than 20. 
if nargout == 1 && numel(P) > 20 && ...
        (exist('vecnorm','builtin') || size(xy,1)==1)
    [~,d]=fminsearch(@(x) vecnorm(xy - repmat([x polyval(P,x)], size(xy,1), 1)), 0, ...
    optimset('tolfun',0));

elseif nargout == 1 && exist('vecnorm','builtin') && size(xy,1) > 10 % TODO find the optimal n in size(xy,1)>n
    [~,d]=fminsearch(@(x) vecnorm(xy - repmat([x polyval(P,x)], size(xy,1), 1)), 0, ...
        optimset('tolfun',0));
    
else
    d = nan(size(xy,1),1);
    p = nan(size(xy));
    for i = 1:size(xy,1)
        xyc = xy(i,:);
        % solve the convolution of (Px-y0) and dPx/dx first
        q = conv([P(1:end-1) P(end)-xyc(2)], polyder(P));
        % then (x-x0)
        q(end-1:end) = q(end-1:end) + [1 -xyc(1)];

        % get roots and calculate y
        x = roots(q); % this is very slow at high degrees.
        y = polyval(P,x);

        % calculate min distance
        [d(i), I] = min(sqrt((x-xyc(1)).^2 + (y-xyc(2)).^2));
        if nargout > 1
            p(i,:) = [x(I),y(I)];
        end
    end
end

% more @
% https://www.mathworks.com/matlabcentral/cody/problems/43642-euclidean-distance-from-a-point-to-a-polynomial