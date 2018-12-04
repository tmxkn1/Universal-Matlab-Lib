function [xc,yc,R,a,me,xn,yn] = circfit(x,y)

[xc,yc,R] = HyperSVD([x,y]);

if nargout >=4
    a = [-2*xc, -2*yc];
    a(3) = (a(1)^2+a(2)^2)/4-R^2;
end

if nargout >= 5
    me = abs( sqrt( (x-xc).^2 + (y-yc).^2 ) - R );
end

if nargout >= 6
    xn = x;
    yn = sign(y-yc).*sqrt( R^2 - (x-xc).^2 ) + yc;
    % for x outside the circle, use y value for calculation
    id = any(imag(yn),2);
    if any(id)
        xn(id) = sign(x(id)-xc).*sqrt( R^2 - (y(id)-yc).^2 ) + xc;
        yn(id) = y(id);
    end
end

function [xc, yc, r, a] = KasaFit(x, y)
%
%   [xc yx R] = circfit(x,y)
%
%   fits a circle  in x,y plane in a more accurate
%   (less prone to ill condition )
%  procedure than circfit2 but using more memory
%  x,y are column vector where (x(i),y(i)) is a measured point
%
%  result is center point (yc,xc) and radius R
%  an optional output is the vector of coeficient a
% describing the circle's equation
%
%   x^2+y^2+a(1)*x+a(2)*y+a(3)=0
%
%  By:  Izhak bucher 25/oct /1991,
x=x(:); y=y(:);
a=[x y ones(size(x))]\-(x.^2+y.^2);
xc = -.5*a(1);
yc = -.5*a(2);
r  =  sqrt((a(1)^2+a(2)^2)/4-a(3));

function [xc, yc, r] = HyperSVD(p)

%--------------------------------------------------------------------------
%  
%     Algebraic circle fit with "hyperaccuracy" (with zero essential bias)
%
%     Input:  XY(n,2) is the array of coordinates of n points x(i)=XY(i,1), y(i)=XY(i,2)
%
%     Output: Par = [a b R] is the fitting circle:
%                           center (a,b) and radius R
%
%     Note: this is a version optimized for stability, not for speed
%
%--------------------------------------------------------------------------

centroid = mean(p);   % the centroid of the data set

X = p(:,1) - centroid(1);  %  centering data
Y = p(:,2) - centroid(2);  %  centering data
Z = X.*X + Y.*Y;
ZXY1 = [Z X Y ones(length(Z),1)];
[~,S,V]=svd(ZXY1,0);
if (S(4,4)/S(1,1) < 1e-12)   %  singular case
    A = V(:,4);
else                         %  regular case
    R = mean(ZXY1);
    N = [8*R(1) 4*R(2) 4*R(3) 2; 4*R(2) 1 0 0; 4*R(3) 0 1 0; 2 0 0 0];
    W = V*S*V';
    [E,D] = eig(W*(N\W));
    [~,ID] = sort(diag(D));
    Astar = E(:,ID(2));
    A = W\Astar;
end

Par = [-(A(2:3))'/A(1)/2+centroid , sqrt(A(2)*A(2)+A(3)*A(3)-4*A(1)*A(4))/abs(A(1))/2];
xc = Par(1);
yc = Par(2);
r = Par(3);
