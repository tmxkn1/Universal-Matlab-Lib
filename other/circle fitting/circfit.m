function   [xc,yc,R,a,me,xn,yn] = circfit(x,y)
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
a=[x y ones(size(x))]\[-(x.^2+y.^2)];
xc = -.5*a(1);
yc = -.5*a(2);
R  =  sqrt((a(1)^2+a(2)^2)/4-a(3));

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
%     
%     if mean( abs(y1-y) ) > mean( abs(y2-y) )
%         y1 = y2;
%     end
%     y = y1;
end