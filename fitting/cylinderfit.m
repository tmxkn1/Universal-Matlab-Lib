function [ax, c, r, eMin, angles] = cylinderfit(p, searchRange, iMax)
if nargin < 2
    searchRange = [-pi/2,pi/2,-pi/2,pi/2];
end

if nargin < 3
    iMax = [(searchRange(2)-searchRange(1))*100,(searchRange(4)-searchRange(3))*100];
    iMax(iMax <= 1) = 1;
end

kMax = iMax(1);
jMax = iMax(2);
eMin = inf;

cmass = mean(p);
p_demass = p - cmass;

omega = linspace(searchRange(1), searchRange(2), kMax);
phi = linspace(searchRange(3), searchRange(4), jMax);
for k = 1:kMax
    for j = 1:jMax
        ax_ = [0,0,1] * rotx(omega(k)) * rotz(phi(j));
        [err, c_, rSqr_] = calculateG(p_demass, ax_);
        if err < eMin
            eMin = err;
            ax = ax_;
            c = c_';
            r = rSqr_;
            angles = [omega(k), phi(j)];
        end
    end
end
c = c + cmass;
r = sqrt(r);
end

function [err, PC, rSqr] = calculateG(X, W)
n = size(X,1);

% projection matrix
P = eye(3) - W'*W;

Y = P * X';

% dot product of each Y
sqrl = sum(Y.^2);

% sum of the outer products of Y
A = zeros(3);
for i = 1:n
    A = A+Y(:,i)*Y(:,i)';
end

% skew matrix
S = [0, -W(3), W(2); W(3), 0, -W(1); -W(2), W(1), 0];
% A hat
Ahat = S*A*S';


PC = Ahat/trace(Ahat*A)*sum(sqrl.*Y,2);
mPC = repmat(PC,1,size(Y,2));

err = mean((sqrl - mean(sqrl) - 2*dot(Y, mPC)).^2);
rSqr = mean(sum((PC - Y).^2));
end