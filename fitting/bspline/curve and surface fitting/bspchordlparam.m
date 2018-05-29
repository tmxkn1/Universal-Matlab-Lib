function [u, cumulength] = bspchordlparam(var1, var2, var3, var4)
% This function applies the chord length method to parameterise a set of
% data points.
% USE:
% returnedData = sub_Bspline_chordL(data, no_of_seg)
%   WHERE
%     data:
%       a n-by-3 matrix; n is the number of total points. Each column of
%       the matrix contains X-, Y- and Z- coordinates of the data points
%       that need to be parameterised.
%
%     no_of_seg:
%       a scalar that defines the number of segments or the maximum value
%       of u.
%
% returnedData = sub_Bspline_chordL(x, y, no_of_seg)
% returnedData = sub_Bspline_chordL(x, y, z, no_of_seg)

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