function setup()
% setup  Add SchizoVis (renamed from SchizoViz) source to MATLAB path.
%   Run this once per session (or add to startup.m) before using examples/tests.

rootDir = fileparts(mfilename('fullpath'));
srcDir  = fullfile(rootDir, 'src');
if exist(srcDir, 'dir')
    addpath(srcDir);
    fprintf('Added src to path: %s\n', srcDir);
else
    warning('src directory not found (%s).', srcDir);
end
end
