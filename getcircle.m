function p = getcircle(c,r,n,v)

if nargin < 4
    v = [1 0 0];
end

[d3, d2] = size(c);
if d2>3
    error('c must be a n-by-3 matrix.')
end
if d3~=numel(r)
    error('The number of hole centres do not equal to the number of radii.');
end

% set z = 0 if 2d
is2d = 0;
if d2==2
    is2d = 1;
    c = [c zeros(size(r))];
    d2 = 3;
end

d1 = round(n/4)*4;
if d1~=n
    warning('n is indivisible by 4, hence rounded to %d', n);
end

p = zeros(d1,d2,d3);

% (x-a)^2 + (y-b)^2 + (z-c)^2 = r^2
func = @(x,y,z) linspace(x,x+y,z);
x = arrayfun(func,c(:,1), r(:), d1/4);
x(end+1:end+d1/4) = fliplr(x);
x(end+1:end+d1/2) = x-r;