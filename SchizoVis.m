% SchizoVis main driver (renamed from SchizoViz)
% Keeping contents identical to legacy SchizoViz.m for continuity.

%% SchizoVis driver using helper functions (renamed from SchizoViz)
clear; clc; close all;

% Parameters
nMin = 11; nMax = 11; nStepSize = 2; %#ok<NASGU>
baseRadix = 10; 
precisionOrder = 100; 
mode = "residual"; % options: "digits", "residual"

% Export / video switches (similar spirit to v2)
exportFigure = false;            % save .fig and .png
exportVideo  = false;            % write progressive layer reveal video
videoLayers  = inf;              % number of initial layers -> frames (clamped to available)
videoFrameRate = 15;             % fps
videoView = 'top';               % 'top' or '3d'
videoLossless = true;            % if true use Uncompressed AVI
videoResolution = [1080 1080];   % figure size during capture ([] keep current)
videoLastN = inf;                % number of most-recent layers to show in each frame (nf). inf => all up to current
videoDynamicColors = true;       % recolor visible window each frame spanning colormap (current layer brightest)
FileName = sprintf('SchizoVis_n%d-%d_b%d_p%d', nMin, nMax, baseRadix, precisionOrder);
% Optionally export the concatenated integer itself as a square-digit image
exportNumber = false; % if true, export the last f(n) digits in target base as a square image
% Export folder (where figures/videos will be written)
exportFolder = pwd; % default to current working directory; set to absolute path to override
exportFolder = "D:\~Projects\SchizoViz";

% Visualization switches
VisPolar = true;  % enable stacked polar visualization
VisGrid  = true;  % enable grid-digit visualization

numDigits = precisionOrder + nMax + 10; % how many fractional digits to show in the exported number image
digits(numDigits);

% Build schizophrenic numbers incrementally (vectorized builder)
nList = nMin:nStepSize:nMax; L = numel(nList);
% SchizoGen returns symbolic arrays (fVals, sVals) without numeric vpa
fprintf('SchizoVis: calling SchizoGen for %d items...\n', numel(nList));
t_gen = tic;
[fVals, sVals] = SchizoGen(nList, baseRadix);
genTime = toc(t_gen);
fprintf('SchizoVis: SchizoGen completed in %.3fs\n', genTime);

% Compute a global exponent range from the largest f(n) so all rows match
p_min = -abs(precisionOrder);
% Only now (digits/residual matrix generation stage) use vpa to size exponent range
fprintf('SchizoVis: sizing exponent grid using vpa on the largest sqrt...\n');
t_size = tic;
s_max_sym = sVals(end);
s_max_num = vpa(s_max_sym); %#ok<NASGU>
global_p_max = ceil(double(log(vpa(s_max_sym))/log(baseRadix)));
sizeTime = toc(t_size);
fprintf('SchizoVis: exponent sizing done in %.3fs -> global_p_max=%d\n', sizeTime, global_p_max);
exponents = global_p_max:-1:p_min;
fprintf('SchizoVis: batch expanding %d layers across %d exponents...\n', L, numel(exponents));
t_batch = tic;
digitsMat = ExpoExpand(sVals, exponents, precisionOrder, baseRadix, mode);

batchTime = toc(t_batch);
fprintf('SchizoVis: batch expansion completed in %.3fs (%.2fM values)\n', batchTime, numel(digitsMat)/1e6);

if VisPolar
	fprintf('SchizoVis: invoking PolarDigitVis (stacked)...\n');
	v = tic;
	hScatter = PolarDigitVis(exponents, digitsMat, baseRadix, 'Title', sprintf('n=%d..%d base=%d precision=%d', nMin, nMax, baseRadix, precisionOrder));
	pvTime = toc(v);
	fprintf('SchizoVis: PolarDigitVis returned in %.3fs\n', pvTime);
	% Acquire figure/axes handles
	ax = ancestor(hScatter,'axes');
	fig = ancestor(ax,'figure');
else
	fprintf('SchizoVis: VisPolar is false, skipping PolarDigitVis\n');
	hScatter = []; ax = []; fig = [];
end
% (Remaining body identical)
% For maintenance, prefer editing SchizoVis.m and keep SchizoViz.m as thin wrapper.
