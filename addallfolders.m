function p = addallfolders(varargin)
%ADDALLFOLDERS add all subfolders in the current directory to MATLAB search 
% path.
%
% addallfolders(F1, F2,...) ignores the specified folders with names F1,
% F2, etc.
%
% addallfolders(M, ...) add folders and subfolders which contains file M. M
% must be specified in fullpath, for example:
%
% M = 'c:/root/file.txt';
%
% p = addallfolders(...) returns all folders that have been added.
%
% Copyright Zhengyi Jiang, 2017.

curpath = what;
fpath = [];

if nargin > 0
    fpath = fileparts(varargin{1});
    % if 
    if ~isempty(fpath)
        cd(fpath);
        varargin{1} = [];
        c = what;
    else
        c = curpath;
    end
else
    c = curpath;
end

p = genpath(c.path);
if ~isempty(varargin)
    p = strsplit(p,';');
    
    for i = 1:numel(varargin)
        if isempty(p)
            return
        end
        p = p(cellfun(@isempty, strfind(p,varargin{i})));
    end
    p = strjoin(p,';');
end
   
% return to orginal path
if ~isempty(fpath)
    cd(curpath.path);
end

% add path
addpath(p)