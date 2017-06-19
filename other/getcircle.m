function [x,y] = getcircle(c,r,n)
%GETCIRCLE calculate points on a given circle
% [X,Y] = getcircle(C,R,N) returns the N numer of X- and Y-coordinates of 
% circles specified by centres C and radii R.
%
% C is an n-by-2 matrix, each row specifies the X- and Y-coordinates of a
% circle.
%
% R is an n-by-1 or 1-by-n array specifies the radii of the circles.
%
% X and Y are N-by-n matrices, each column of which represents the X- or
% Y-coordinates of a circle.
%
% by Zhengyi Jiang, ELE Advanced Techonologies Ltd., 2017

[d1, d2] = size(c);
if d2~=2
    error('c must be a n-by-2 matrix.')
end
if d1~=numel(r)
    error('The number of hole centres do not equal to the number of radii.');
end

th = linspace(0,2*pi,n+1)';
x = r.*cos(th)+c(:,1)'; 
y = r.*sin(th)+c(:,2)';