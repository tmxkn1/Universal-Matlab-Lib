function val = transformby(val, transmat)

n = sqrt(numel(transmat));
if numel(val(1,:)) ~= n
    val = [val ones(size(val,1),1)];
    val = val * transmat;
    val = val(:,1:n-1);
else
    val = val * transmat;
end

