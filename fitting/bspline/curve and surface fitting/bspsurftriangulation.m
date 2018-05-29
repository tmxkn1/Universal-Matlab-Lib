function [u,w] = bspsurftriangulation(sectionDataInCell)

pc = sectionDataInCell;
noc = numel(pc);

% w_ = 0:1/(noc-1):1;
r=cell2mat(transpose(cellfun(@(x) x(1,:),sectionDataInCell,'UniformOutput', false)));
w_ = bspchordlparam(r,1);

for i = 1:noc
    uc{i} = bspchordlparam(pc{i},1);
    wc{i} = w_(i)*ones(size(uc{i}));
end

u = cell2mat(uc);
w = cell2mat(wc);
u = u(:);
w = w(:);