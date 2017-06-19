function codepad

a=magic(6000);
a=reshape(a,12000000,3);
b=[1 21 42];


%c = a.^b;
dist=sqrt((a(:,1)-b(1)).^2+(a(:,2)-b(2)).^2+(a(:,3)-b(3)).^2);
%dist = sum(bsxfun(@minus, b, a) .^ 2, 2);