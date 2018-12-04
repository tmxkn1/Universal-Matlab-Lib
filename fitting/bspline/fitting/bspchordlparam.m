function [u, cumulength] = bspchordlparam(var1, var2, var3, var4)
% Parametrises a 2/3-D curve using the chord length method.
%
% -------------------------------------------------------------------------
% USE:
%
% O = bspchordlparam(D, S) calculates a set of parametrised coordinates,
%       O, for a 2/3-D curve defined by an n*3/n*2 matrix, D. S is
%       is the number of segments, which is indicates the maximum value of 
%       O.
%
% O = bspchordlparam(X, Y, S) calculates a set of parametrised coordinates,
%       O, for a 2-D curve defined by a set of X- and Y-coordinates.
%
% O = bspchordlparam(X, Y, Z, S) calculates a set of parametrised
%       coordinates, O, for a 3-D curve defined by a set of X-, Y- and 
%       Z-coordinates.
%
% [O, C] = bspchordlparam(...) also returns the cumulative length of the 
%       curve.
%
% -------------------------------------------------------------------------
% See also: bspbasis, bspcontrolpt, bspgenpoints, bspknot, bspcurvefit

%% Validate input
if nargin == 2
    data = var1;
    no_of_seg = var2;
    [a, b] = size(data);
    if b~=3
        error('Input data should be a n-by-3 matrix.');
    elseif a<2
        error('Insufficient data.');
    end
    if ~isscalar(no_of_seg)
        error('Parameter no_of_seg must be a scalar.');
    end    
    x = data(:,1); y = data(:,2); z = data(:,3);
elseif nargin == 3
    x = var1;
    y = var2;
    no_of_seg = var3;
    z = zeros(size(x));
elseif nargin == 4
    x = var1;
    y = var2;
    z = var3;
    no_of_seg = var4;
end

%% The chord length method
% Sort out data for calculation
np = size(data,1);

% Calculate all the small chord lengths.
uvalue = [0; sqrt(diff(x).^2 +diff(y).^2+diff(z).^2)];
cumulength = cumsum(uvalue);
u = transpose(cumulength/cumulength(np)*no_of_seg);