startupparam.ignoreddir = {'.git'};

startupparam.f = mfilename('fullpath');
[startupparam.path_, ~, ~] = fileparts(startupparam.f);
startupparam.p = strsplit(genpath(startupparam.path_),';');

for i = 1:numel(startupparam.ignoreddir)
    if isempty(startupparam.p)
        break;
    end
    startupparam.p = startupparam.p(cellfun(@isempty, strfind(startupparam.p,startupparam.ignoreddir{i})));
end

if ~isempty(startupparam.p)
    addpath(strjoin(startupparam.p,';'))
end

clear startupparam i;