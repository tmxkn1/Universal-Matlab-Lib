function [CP,CLU,CRU,CLW,CRW] = bezsubdiv(C,U,W)
%BEZSUBDIV Subdivision of a Bezier curve or a Bezier patch.
%
% CP = BEZSUBDIV(C,U) splits a Bezier curve with control points C at 
% parameter U, and returns new control points CP.
% - C is a k-by-3 matrix, where k is the order of the curve
% - CP is a kp-by-3 matrix, where kp = 2*k-1
%
% [CP, CL, CR] = BEZSUBDIV(...) also returns k-by-3 matrices CL and CR for
% control points CL and CR for left and right segments.
%
% CP = BEZSUBDIV(C,U,W) splits a Bezier patch with control points C at
% parameter U and W, and returns new control points CP.
% - C is a ku-by-kw-3 matrix, where ku is the order of the patch in the
% u-direction and kw is that in the w-direction.
% - CP is a kup-by-kwp-3 matrix, where kup = 2*ku-1 and kwp = 2*kw-1
% 
% [CP,CLU,CRU,CLW,CRW] = BEZSUBDIV(...) also returns ku-by-kw-3 matrices 
% CLU, CRU, CLW and CRW for control points for top left, top right, bottom 
% left and bottom right patches.

if nargin == 2 % Bezier curve
    CLU = C(1,:); 
    CRU = C(end,:);
    for i = 2:size(C,1)
       C = (1-U)*C(1:end-1,:) + U*C(2:end,:);
       CLU = [CLU; C(1,:)]; 
       CRU = [C(end,:); CRU];
    end
    CLW = 0; CRW = 0;
    CP = [CLU; CRU(2:end,:)];
    
else % Bezier surface
    % split in u-direction first.
    CLU = C(1,:,:); 
    CRU = C(end,:,:);
    for i = 2:size(C,1)
       C = (1-U)*C(1:end-1,:,:) + U*C(2:end,:,:);
       CLU = [CLU; C(1,:,:)]; 
       CRU = [C(end,:,:); CRU];
    end

    % then split in w-direction
    % split clu
    C = CLU;
    CLU = C(:,1,:); 
    CLW = C(:,end,:);
    for i = 2:size(C,1)
       C = (1-W)*C(:,1:end-1,:) + W*C(:,2:end,:);
       CLU = [CLU  C(:,1,:)]; 
       CLW = [C(:,end,:)  CLW];
    end
    
    % split cru
    C = CRU;
    CRU = C(:,1,:); 
    CRW = C(:,end,:);
    for i = 2:size(C,1)
       C = (1-W)*C(:,1:end-1,:) + W*C(:,2:end,:);
       CRU = [CRU  C(:,1,:)]; 
       CRW = [C(:,end,:)  CRW];
    end
    
    % combine control points
    CP = [CLU CLW(:,2:end,:); CRU(2:end,:,:) CRW(2:end,2:end,:)];
end