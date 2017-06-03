function [CN, TN] = bsp2bez(var1, var2, var3, var4, var5)
%BSP2BEZ Decomposes a clamped B-spline curve or surface into several Bezier
% curves or patches.
%
% CN = BSP2BEZ(T, C, K) decomposes a K-th order clamped B-spline curve with
% knots T and control points C into Bezier curves with control points CN.
% - C is a N-by-3 matrix, where N by law equals numel(T)-K.
% 
% CN = BSP2BEZ(TU, TW, C, KU, KW) decomposes a clamped B-spline curve of
% order KU with knots TU in the u-direciton and of order KW with knots TW 
% in the w-direction into Bezier patches with control points CN.
% - C is a NU-by-NW-by-3 matrix, where NU = numel(TU)-KU and 
%   NW = numel(TW)-KW.
%
% [CN, TN] = bsp2bez(...) also returns the knots after insertion.
%
% ref: http://www.infogoaround.org/JBook/bstobez.html

if nargin == 3 % curve conversion
    ts = var1(:);
    c = var2;
    k = var3;
    
    % give initial values to new knots and control points
    TN = ts;
    d = c;
    r = k;
    while r <numel(TN)-k
        % make a copy
        d_ = d;
        % update the number of control points
        m = numel(d(:,1));
        t = TN(r+1);
        s = sum(TN == t);
        % number of knots to insert
        p = k-1-s;
        if p < 0 
            continue;
        end

        % compute k+p-2 affected control points
        b(1:k,1,1) = d(r-k+2:r+1,1);
        b(1:k,1,2) = d(r-k+2:r+1,2);
        b(1:k,1,3) = d(r-k+2:r+1,3);
        for j = 1:p
            % code in reference:
            %{
            for i = j:k-1
                alpha = (t - tsnew(r-k+i+2)) / (tsnew(r+2+i-j) - tsnew(r-k+i+2));
                b(i+1,j+1,:) = b(i,j,:)+alpha*(b(i+1,j,:)-b(i,j,:));
            end
            %}
            % to matrix form:
            alpha = (t - TN(r-k+2+j:r+1)) ./ (TN(r+2:r-j+k+1) - TN(r-k+j+2:r+1));
            b(j+1:k,j+1,:) = b(j:k-1,j,:)+alpha.*(b(j+1:k,j,:)-b(j:k-1,j,:));
        end

        % copy the diagonal points
        di = [diag(b(:,:,1)) diag(b(:,:,2)) diag(b(:,:,3))];
        d(r-k+3:r-k+2+p,:) = di(2:end,:); 

        % copy the pth row points
        d(r-k+3+p:r-k+2+2*p,:) = squeeze(b(p+1, p:-1:1, :)); 

         % copy final part of unaffected ponts
        d(2*p-k+r-s+4:2*p-k+m+2,:) = d_(r-s+2:m,:);

        % insert knots
        TN = [TN(1:r+1); t*ones(p,1); TN(r+2:end)];

        % plot(d(r-k+3:r-k+2+2*p,1),d(r-k+3:r-k+2+2*p,2),'o-');
        % increment r
        r = r+k-1;
    end
    CN = d;
elseif nargin == 5 % surface conversion
    % give initial values to new knots and control points
    ts = {var1(:), var2(:)};
    kuw = [var4,var5];
    CN = var3;
    % go through the control points along u and w direction respectively
    % dir == 1 means u
    % dir == 2 means w
    for dir = 1:2
        % order in the current direction
        k = kuw(dir);
        for id = 1:size(CN,3-dir) 
            if dir == 1
                d = squeeze(CN(:,id,:));
            else
                d = squeeze(CN(id,:,:));
            end
            tstemp = ts {dir};
            r = k;
            
            % the while loop is identical to previous.
            while r <numel(tstemp)-k
                % make a copy
                d_ = d;
                % update the number of control points
                m = numel(d(:,1));
                t = tstemp(r+1);
                s = sum(tstemp == t);
                % number of knots to insert
                p = k-1-s;
                if p < 0 
                    continue;
                end

                % compute k+p-2 affected control points
                b(1:k,1,1) = d(r-k+2:r+1,1);
                b(1:k,1,2) = d(r-k+2:r+1,2);
                b(1:k,1,3) = d(r-k+2:r+1,3);
                for j = 1:p
                    alpha = (t - tstemp(r-k+2+j:r+1)) ./ (tstemp(r+2:r-j+k+1) - tstemp(r-k+j+2:r+1));
                    b(j+1:k,j+1,:) = b(j:k-1,j,:)+alpha.*(b(j+1:k,j,:)-b(j:k-1,j,:));
                end

                % copy the diagonal points
                di = [diag(b(:,:,1)) diag(b(:,:,2)) diag(b(:,:,3))];
                d(r-k+3:r-k+2+p,:) = di(2:end,:); 

                % copy the pth row points
                d(r-k+3+p:r-k+2+2*p,:) = squeeze(b(p+1, p:-1:1, :)); 

                 % copy final part of unaffected ponts
                d(2*p-k+r-s+4:2*p-k+m+2,:) = d_(r-s+2:m,:);

                % insert knots
                tstemp = [tstemp(1:r+1); t*ones(p,1); tstemp(r+2:end)];
                
                % plot3(d(r-k+3:r-k+2+2*p,1),d(r-k+3:r-k+2+2*p,2),d(r-k+3:r-k+2+2*p,3),'o-k');
                % increment r
                r = r+k-1;
            end
            cnewtemp{id} = reshape(d, size(d,1),1,3);
        end
        % after control points are calculated for one of the directions
        CN = cell2mat(cnewtemp);
        TN{dir} = tstemp;
        % clear b
        b=[];
    end
end