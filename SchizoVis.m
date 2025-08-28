%% SchizoVis driver using helper functions
clear; clc; close all;

% Parameters
nMin = 001; nMax = 101; nStepSize = 2; %#ok<NASGU>
baseRadix = 10; 
precisionOrder = 500; 
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
% Export folder (where figures/videos will be written)
exportFolder = pwd; % default to current working directory; set to absolute path to override
exportFolder = "D:\~Projects\SchizoVis";

% Visualization switches
VisPolar = true; % show polar stacked visualization
VisGrid  = false;  % show grid digit visualization (last sqrt only for static; sweep video if exportVideo)

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
	tv = tic;
	hScatter = PolarDigitVis(exponents, digitsMat, baseRadix, 'Title', sprintf('n=%d..%d base=%d precision=%d', nMin, nMax, baseRadix, precisionOrder));
	pvTime = toc(tv);
	fprintf('SchizoVis: PolarDigitVis returned in %.3fs\n', pvTime);
	% Acquire figure/axes handles
	ax = ancestor(hScatter,'axes');
	fig = ancestor(ax,'figure');
else
	hScatter = []; ax = []; fig = [];
end

% Grid visualization (non-export) when VisGrid is enabled and exportNumber not used
if VisGrid && ~exportVideo
	try
		GridDigitVis(sVals, baseRadix, '', nList, 'FractionDigits', precisionOrder, 'Mode', mode, 'Transparent', true, 'Show', true);
		try
			figGrid = gcf; set(figGrid,'WindowStyle','docked');
		catch
		end
	catch ME
		warning(ME.identifier, 'Could not create multi-row grid visualization: %s', ME.message);
	end
elseif VisGrid && exportVideo
fprintf('SchizoVis: exporting grid video (one frame per n) using GridDigitVis...\n');
	if ~isfolder(exportFolder); mkdir(exportFolder); end
	gridVideoName = sprintf('%s_grid', FileName);
	proceed = true; extG = '.avi';
	try
		if videoLossless
			vwGrid = VideoWriter(fullfile(exportFolder, [gridVideoName extG]), 'Uncompressed AVI');
		else
			try
				extG = '.mp4'; vwGrid = VideoWriter(fullfile(exportFolder, [gridVideoName extG]), 'MPEG-4');
			catch
				extG = '.avi'; vwGrid = VideoWriter(fullfile(exportFolder, [gridVideoName extG]), 'Motion JPEG AVI');
			end
		end
		vwGrid.FrameRate = min(videoFrameRate,30);
		open(vwGrid);
	catch ME
		warning(ME.identifier,'Grid video: could not open writer: %s', ME.message);
		proceed = false;
	end
	if proceed
		for idxN = 1:numel(nList)
			try
				GridDigitVis(sVals(idxN), baseRadix, '', nList(idxN), 'FractionDigits', precisionOrder, 'Mode', mode, 'Transparent', true, 'Show', true);
				figGV = gcf;
				% black background for consistency
				try
					set(figGV,'Color','k');
				catch
				end
				frameGV = getframe(figGV);
				writeVideo(vwGrid, frameGV);
				close(figGV);
				if idxN==1 || idxN==numel(nList) || mod(idxN, max(1,round(numel(nList)/10)))==0
					fprintf('  Grid frame %d/%d (%.1f%%)\n', idxN, numel(nList), 100*idxN/numel(nList));
				end
			catch MEf
				warning(MEf.identifier,'Grid video: frame %d failed: %s', idxN, MEf.message);
			end
		end
		try
			close(vwGrid);
		catch
		end
		fprintf('SchizoVis: grid video export complete (%s).\n', [gridVideoName extG]);
	end
end

% -------- Figure export --------
if exportFigure
	try
		if ~isfolder(exportFolder); mkdir(exportFolder); end

		% Export existing polar figure only (do not recreate)
		if VisPolar
			if exist('fig','var') && isgraphics(fig,'figure')
				figPathFig = fullfile(exportFolder, sprintf('%s_polar.fig', FileName));
				figPathPng = fullfile(exportFolder, sprintf('%s_polar.png', FileName));
				savefig(fig, figPathFig);
				exportgraphics(fig, figPathPng, 'Resolution', 300);
				fprintf('SchizoVis: polar figure exported -> %s , %s\n', figPathFig, figPathPng);
			else
				warning('SchizoVis:NoPolarFig','Polar figure handle not found; skipping polar export.');
			end
		end

		% Export existing grid figure only (do not recreate)
		if VisGrid
			% prefer figGrid then figGV if any exist
			gridFigHandle = [];
			if exist('figGrid','var') && isgraphics(figGrid,'figure')
				gridFigHandle = figGrid;
			elseif exist('figGV','var') && isgraphics(figGV,'figure')
				gridFigHandle = figGV;
			end

			if ~isempty(gridFigHandle)
				gridPathFig = fullfile(exportFolder, sprintf('%s_grid.fig', FileName));
				gridPathPng = fullfile(exportFolder, sprintf('%s_grid.png', FileName));
				try
					savefig(gridFigHandle, gridPathFig);
					exportgraphics(gridFigHandle, gridPathPng, 'Resolution', 150);
					fprintf('SchizoVis: grid figure exported -> %s , %s\n', gridPathFig, gridPathPng);
				catch ME2
					warning(ME2.identifier,'Grid figure export failed: %s', ME2.message);
				end
			else
				warning('SchizoVis:NoGridFig','Grid figure handle not found; skipping grid export.');
			end
		end

	catch ME
		warning(ME.identifier,'Figure export failed: %s', ME.message);
	end
end

% -------- Video export (progressive reveal) --------
if exportVideo && VisPolar && ~isempty(hScatter)
	vLayers = min(videoLayers, L);
	fprintf('SchizoVis: exporting video (%d layers -> frames) to %s ...\n', vLayers, [FileName '_polar']);
	% VideoWriter setup
	if videoLossless
		ext = '.avi'; vw = VideoWriter(fullfile(exportFolder, [FileName '_polar' ext]), 'Uncompressed AVI');
	else
		try
			ext = '.mp4';
			vw = VideoWriter(fullfile(exportFolder, [FileName '_polar' ext]), 'MPEG-4');
		catch
			ext = '.avi';
			vw = VideoWriter(fullfile(exportFolder, [FileName '_polar' ext]), 'Motion JPEG AVI');
		end
	end
	vw.FrameRate = videoFrameRate; open(vw);
	origPos = get(fig,'Position'); origFigColor = get(fig,'Color'); resized=false;
	% use black background for video frames
	try
		set(fig,'Color','k');
		set(ax,'Color','k');
	catch
	end
	if ~isempty(videoResolution) && numel(videoResolution)==2
		try
			newPos = origPos; newPos(3:4) = videoResolution;
			set(fig,'Position',newPos);
			drawnow;
			resized=true;
		catch ME
			warning(ME.identifier,'Could not resize figure: %s', ME.message);
		end
	end

	% Derive clockR & zRange similar to PolarDigitVis for consistent axis locking (only if polar exists)
	gmax = exponents(1);
	if gmax > 0
		r_vals = 1 - exponents / gmax;
	else
		emin = exponents(end);
		if emin == 0
			r_vals = ones(size(exponents));
		else
			r_vals = (exponents - emin)/(0 - emin);
		end
	end
	clockR = max(r_vals);
	L = size(digitsMat,1); zRange = [0 clockR];
	padFrac = 0.03; xr = 2*clockR; pad = padFrac * xr;
	% Only set axis limits and extract scatter data if Polar visualization was created
	if ~isempty(ax) && isgraphics(ax,'axes') && ~isempty(hScatter)
		if clockR>0
			xlim(ax, [-clockR - pad, clockR + pad]);
			ylim(ax, [-clockR - pad, clockR + pad]);
		else
			xlim(ax,[-1 1]); ylim(ax,[-1 1]);
		end
		zpad = max(0.05*clockR, 0.05);
		zlim(ax, [zRange(1)-zpad, zRange(2)+zpad]);
		try
			ax.XLimMode='manual';
			ax.YLimMode='manual';
			ax.ZLimMode='manual';
		catch
			axis(ax,'manual');
		end
		drawnow;
		% Flatten data for later video masking
		allX = get(hScatter,'XData');
		allY = get(hScatter,'YData');
		allZ = get(hScatter,'ZData');
		M = numel(exponents); % per-layer column count
	else
		% No polar scatter: create safe placeholders so later code can check emptiness
		allX = []; allY = []; allZ = [];
		M = numel(exponents);
	end

	% Flatten data for later video masking
	allX = get(hScatter,'XData');
	allY = get(hScatter,'YData');
	allZ = get(hScatter,'ZData');
	M = numel(exponents);

	% Prepare masked arrays (start empty)
	xMask = nan(size(allX)); yMask = nan(size(allY)); zMask = nan(size(allZ));
	if videoDynamicColors
		cmapVid = colormap(ax);
		cRows = size(cmapVid,1);
	end
	for li = 1:vLayers
		% sliding window: show layers [startLayer .. endLayer] where endLayer=li
		endLayer = min(li, L);
		startLayer = max(1, li - min(videoLastN, L) + 1);
		% compute flattened index range to show
		rngShow = (startLayer-1)*M + 1 : endLayer*M;
		% clear mask then fill only the window
		xMask(:) = NaN; yMask(:) = NaN; zMask(:) = NaN;
		xMask(rngShow) = allX(rngShow); yMask(rngShow) = allY(rngShow); zMask(rngShow) = allZ(rngShow);
		if videoDynamicColors
			visibleLayers = startLayer:endLayer;
			w = numel(visibleLayers);
			% map oldest -> lowest colormap index, current -> highest index across full colormap
			if w == 1
				layerColorIdx = cRows;
			else
				layerColorIdx = round(linspace(1, cRows, w));
			end
			frameC = zeros(numel(allX),3);
			for kLayer = 1:w
				layer = visibleLayers(kLayer);
				idxRange = (layer-1)*M + 1 : layer*M;
				frameC(idxRange,:) = repmat(cmapVid(layerColorIdx(kLayer),:), M,1);
			end
			% non-visible points get a dark base color
			maskZero = all(frameC==0,2);
			if any(maskZero)
				frameC(maskZero,:) = 0.12; % dark gray
			end
			set(hScatter,'CData',frameC);
		end
		set(hScatter,'XData',xMask,'YData',yMask,'ZData',zMask);
		switch lower(videoView)
			case 'top'; view(ax,2);
			otherwise; view(ax,3);
		end
		title(ax, sprintf('n=%d base=%d precision=%d', nList(li), baseRadix, precisionOrder),'Color','w');
		drawnow limitrate;
		writeVideo(vw, getframe(fig));
		if li==1 || li==vLayers || mod(li, max(1,round(vLayers/10)))==0
			fprintf('  Frame %d/%d (%.1f%%)\n', li, vLayers, 100*li/vLayers);
		end
	end
	close(vw);
	if resized
		try
			set(fig,'Position',origPos);
			drawnow;
		catch ME
			warning(ME.identifier,'Could not restore original figure position: %s', ME.message);
		end
	end
	fprintf('SchizoVis: video export complete (%s).\n', [FileName '_polar' ext]);
end

% (Optionally could restore full scatter data if modified during video) 
if VisPolar && ~isempty(hScatter)
	try
		set(hScatter,'XData',allX,'YData',allY,'ZData',allZ);
	catch ME
		warning(ME.identifier,'Could not restore scatter data: %s', ME.message);
	end
end


