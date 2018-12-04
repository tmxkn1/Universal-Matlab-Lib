function [X,Y,Z] = bezfit(var1, var2, var3, var4, var5)
%BEZFIT Fit Bezier curve(s) or patch(es).
%
% [X,Y,Z] = BEZFIT(C, K) calculates K-th order bezier curve(s) with control
% points C and returns X, Y and Z-coordinates of the fitted bezier curve(s)
% calculated with 101 uniform parameters.
% C is a N-by-3 matrix, where N is the number of segments multiplied by
% K-1.
% 
% [X,Y,Z] = BEZFIT(C, K, U) calculates K-th order bezier curve(s) with 
% control points C and returns X, Y and Z-coordinates of the fitted bezier 
% curve(s) calculated with parameter U.
% 
% [X,Y,Z] = BEZFIT(C, KU, KW) calculates bezier patch(es) of order KU in
% u-direction and order KW in w-direction with control points C and returns
% X, Y and Z-coordinates of the fitted patches with 101 uniform parameters 
% in both direction.
% - C is a NU-by-NW-by-3 matrix, where NU is the the number of patches in 
%  u-direction multiplied by KU-1 and NW is that in w-direction.
%
% [X,Y,Z] = BEZFIT(..., U) uses U for parameters in the u direciton to fit
% the Bezier curve(s) or patch(es)
%
% [X,Y,Z] = BEZFIT(..., U, W) uses U and W for parameters in both
% directions to fit the Bezier curve(s) or patch(es).


ctrlPt = var1;
if size(ctrlPt,3) == 1 % bezier curve
    k = var2;
    % generate uniform u values if not given as an input
    if nargin == 2
        var3 = 0:0.01:1;
    end
    uarray = var3;
    
    % number of u values
    nu = numel(uarray);
    
    % number of segments
    nos = size(ctrlPt,1);
    nos = (nos-1)/(k-1);
    
    % construct matrices
    u = repmat(uarray(:),1,k);
    px = repmat(ctrlPt(:,1)', nu, 1);
    py = repmat(ctrlPt(:,2)', nu, 1);
    pz = repmat(ctrlPt(:,3)', nu, 1);
    
    % calculate bernstein polynomials
    fu = bernsteinPoly(k,u);
    
    % fit curves
    for i = 1:nos
        st = (i-1)*(k-1)+1; ed = st+k-1;
        stx = (i-1)*nu+1; edx = stx+nu-1;
        X(stx:edx,:) = fu*ctrlPt(st:ed,1);
        Y(stx:edx,:) = fu*ctrlPt(st:ed,2);
        Z(stx:edx,:) = fu*ctrlPt(st:ed,3);
    end
    
else % bezier patch
    ku = var2;
    kw = var3;
    % generate uniform u and w values if not given as inputs
    if nargin == 3
        var4 = 0:0.01:1;
        var5 = var4;
    elseif nargin == 4
        var5 = 0:0.01:1;
    end
    uarray = var4;
    warray = var5;
    
    % number of u and w values
    nu = numel(uarray);
    nw = numel(warray);
    
    % construct matrices
    u = repmat(uarray(:),1,ku);
    w = repmat(warray(:),1,kw);
    
    % number of patches in both direction
    [nopu, nopw, ~] = size(ctrlPt);
    nopu = (nopu-1)/(ku-1);
    nopw = (nopw-1)/(kw-1);
    
    % calculate bernstein polynomials
    fu = bernsteinPoly(ku,u);
    fw = bernsteinPoly(kw,w);
    
    % fit patches
    for i = 1:nopw
        for j = 1:nopu
            stu = (j-1)*(ku-1)+1; edu = stu+ku-1;
            stw = (i-1)*(kw-1)+1; edw = stw+kw-1;
            stxu = (j-1)*nu+1; edxu = stxu+nu-1;
            stxw = (i-1)*nw+1; edxw = stxw+nw-1;
            X(stxu:edxu,stxw:edxw) = fu*ctrlPt(stu:edu,stw:edw,1)*fw';
            Y(stxu:edxu,stxw:edxw) = fu*ctrlPt(stu:edu,stw:edw,2)*fw';
            Z(stxu:edxu,stxw:edxw) = fu*ctrlPt(stu:edu,stw:edw,3)*fw';
        end
    end
end

function b = bernsteinPoly(k,param)
b = factorial(k-1)./(factorial(0:k-1).*factorial(k-1:-1:0)).*(1-param).^(k-1:-1:0).*param.^(0:k-1);