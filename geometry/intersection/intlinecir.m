function [x,y] = intlinecir(l,c,r)
%INTLINECIR intersection of multiple lines and a circle
%
% [X,Y] = intlinecir(L,C,R) returns X- and Y-coordinates of intersection
% points between lines defined by y=L(:,1)*x+L(:,2) and a circle whose
% centre is C and radius is R.
%
% by Zhengyi Jiang, ELE Advanced Technologies Ltd, 2017


d = l(:,2)-c(2);

% coefficients of 2nd order equation
a_ = l(:,1).^2 + 1;
b_ = 2*l(:,1).*d - 2*c(1);
c_ = d.^2 + c(1)^2 - r^2;

% solve the equation
d_ = sqrt(b_.^2-4*a_.*c_);

x1 = (-b_ + d_) / 2 ./ a_;
y1 = l(:,1) .* x1 + l(:,2);
x2 = (-b_ - d_) / 2 ./ a_;
y2 = l(:,1) .* x2 + l(:,2);

x = [x1(:) x2(:)]; y = [y1(:) y2(:)];