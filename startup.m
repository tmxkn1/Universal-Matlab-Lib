ignoreddir = {'.git'};

f = mfilename('fullpath');
[path_, ~, ~] = fileparts(f);
p = strsplit(genpath(path_),';');

for i = 1:numel(ignoreddir)
    if isempty(p)
        break;
    end
    p = p(cellfun(@isempty, strfind(p,ignoreddir{i})));
end

if ~isempty(p)
    addpath(strjoin(p,';'))
end

clear;