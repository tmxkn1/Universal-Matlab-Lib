f = mfilename('fullpath');
[path_, ~, ~] = fileparts(f);
addpath(genpath(path_))

clear;