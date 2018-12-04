function [u,w] = bspsurfchordlparam(sectionDataInCell)
% Parametrises a set of surface points.
%
% -------------------------------------------------------------------------
% USE:
%
% [U, W] = bspsurfchordlparam(S) parametrises the surface points, S, using
% the chord-length method, and returns U and W, which are parametrised 
% coordinates in the u- and w-direction respectively.
%
% S is a cell, each element of the cell contains the X-, Y- and 
% Z-coordinates of one section of the surface.
%
% -------------------------------------------------------------------------
% See also: bspsurfspparam, bspsurffit


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