function setup()
%SETUP Adds the project's source directory to the MATLAB path.
%   SETUP should be run once per MATLAB session before using any of the
%   project's functions or running the examples. It locates the 'src'
%   directory relative to its own path and adds it to the top of the
%   MATLAB search path.
%
%   This ensures that all the core functions of the SchizoVis project,
%   such as SchizoGen and PolarDigitVis, are available to be called from
%   the MATLAB command window or from your own scripts.
%
%   Usage:
%       From the root directory of the project in MATLAB, simply run:
%       >> setup

rootDir = fileparts(mfilename('fullpath'));
srcDir  = fullfile(rootDir, 'src');
if exist(srcDir, 'dir')
    addpath(srcDir);
    fprintf('Added src to path: %s\n', srcDir);
else
    warning('src directory not found (%s).', srcDir);
end
end
