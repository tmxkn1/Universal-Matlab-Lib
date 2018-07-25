function val = transformby(val, transmat, doTranslate)
%

if nargin < 3
    doTranslate = 1;
end

n = sqrt(numel(transmat));
if numel(val(1,:)) ~= n
    val = [val ones(size(val,1),1)];
    
    if ~doTranslate
        transmat(end,1:numel(val(1,:))-1) = 0;
        transmat(end) = 1;
    end
    
    val = val * transmat;
    val = val(:,1:n-1);
else
    val = val * transmat;
end

