%% SchizoVis driver using helper functions
clear; clc; close all;

% Parameters
nMin = 001; nMax = 101; nStepSize = 2; %#ok<NASGU>
baseRadix = 2; 
precisionOrder = 01000; 
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
VisPolar = false; % show polar stacked visualization
VisGrid  = true;  % show grid digit visualization (last sqrt only for static; sweep video if exportVideo)

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
	fprintf('SchizoVis: VisPolar is false, skipping PolarDigitVis\n');
	hScatter = []; ax = []; fig = [];
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

% (Removed legacy exportNumber feature to reduce switch conflicts)

% -------- Grid digit sweep video (schizo number across n) --------
% Control grid video generation with the top-level exportVideo flag so both
% polar and grid videos obey the same setting. Set exportVideo=false to skip all videos.
% Grid sweep video only if both exportVideo and VisGrid and multiple n values
if exportVideo && VisGrid && numel(nList) > 1
	fprintf('SchizoVis: exporting grid digit sweep video across n values...\n');
	gridFractionDigits = min(precisionOrder, 500); % cap for performance
	gridVideoName = sprintf('%s_gridSweep', FileName);
	if ~isfolder(exportFolder); mkdir(exportFolder); end
	% Determine max integer digit count using last schizo number approx
	lastNumApprox = vpa(sVals(end), 40);
	if lastNumApprox > 0
		maxIntDigits = max(1, floor(double(log(lastNumApprox)/log(baseRadix))) + 1);
	else
		maxIntDigits = 1;
	end
	maxTotalDigits = maxIntDigits + gridFractionDigits;
	side = ceil(sqrt(maxTotalDigits)); rowsG = side; colsG = side; frameDigitsCapacity = rowsG*colsG;
	% Video writer
	try
		if videoLossless
			vwGrid = VideoWriter(fullfile(exportFolder, [gridVideoName '.avi']), 'Uncompressed AVI');
		else
			try
				vwGrid = VideoWriter(fullfile(exportFolder, [gridVideoName '.mp4']), 'MPEG-4');
			catch
				vwGrid = VideoWriter(fullfile(exportFolder, [gridVideoName '.avi']), 'Motion JPEG AVI');
			end
		end
		vwGrid.FrameRate = min(videoFrameRate, 30);
		open(vwGrid);
	catch ME
		warning(ME.identifier,'Could not open grid sweep video writer: %s', ME.message);
		exportGridVideo = false;
	end
	if exportGridVideo
	% Figure for grid video (black background for video).
	% Use docked window for grid visualization to keep it attached to MATLAB desktop.
	figG = figure('Visible','off','Color','k','Units','pixels','Position',[120 120 800 800], 'WindowStyle', 'docked');
	axG = axes('Parent',figG); axis(axG,'off'); hold(axG,'on');
	set(axG,'Color','k');
	cmapG = parula(baseRadix);
	colormap(axG, cmapG);
	caxis(axG,[0 baseRadix-1]);
	cbG = colorbar(axG,'southoutside'); cbG.Ticks = 0:baseRadix-1;
	lbl = strings(baseRadix,1); for di=0:baseRadix-1; if di<10, lbl(di+1)=string(di); else, lbl(di+1)=char('A'+(di-10)); end; end; cbG.TickLabels=cellstr(lbl);
	cbG.Color = [1 1 1]; cbG.Label.Color = [1 1 1];
	hImg = imagesc(axG, nan(rowsG, colsG)); set(axG,'YDir','reverse'); axis(axG,'image');
		for idxN = 1:numel(nList)
			symVal = sVals(idxN);
			% integer digit count approx
			intPart = floor(symVal);
			intDigits = [];
			if intPart == 0
				intDigits = 0;
			else
				tmp = intPart; loopCap=1e5; it=0; baseSym = sym(baseRadix);
				while tmp > 0 && it < loopCap
					[q,r] = quorem(tmp, baseSym); intDigits = [double(r) intDigits]; %#ok<AGROW>
					tmp = q; it=it+1;
				end
			end
			% fractional digits
			fracPart = symVal - floor(symVal); fracDigits = [];
			if gridFractionDigits>0 && fracPart~=0
				decPrec = ceil(gridFractionDigits * log(baseRadix)/log(10)) + 8;
				fnum = vpa(fracPart, decPrec);
				for kfd=1:gridFractionDigits
					fnum = fnum * baseRadix; d = floor(fnum); fracDigits(end+1)=double(d); %#ok<AGROW>
					fnum = fnum - d;
				end
			end
			allDigitsFrame = [intDigits fracDigits];
			if numel(allDigitsFrame) > frameDigitsCapacity
				allDigitsFrame = allDigitsFrame(1:frameDigitsCapacity); % truncate
			end
			padNeeded = frameDigitsCapacity - numel(allDigitsFrame);
			if padNeeded>0
				allDigitsFrame = [allDigitsFrame nan(1,padNeeded)];
			end
			imgFrame = reshape(allDigitsFrame, colsG, rowsG)';
			% update image
			hImg.CData = imgFrame;
			hImg.AlphaData = ~isnan(imgFrame);
			title(axG, sprintf('Schizo sqrt(f(%d)) base %d frac %d', nList(idxN), baseRadix, gridFractionDigits),'Interpreter','none','Color','w');
			drawnow limitrate;
			writeVideo(vwGrid, getframe(figG));
			if idxN==1 || idxN==numel(nList) || mod(idxN, max(1,round(numel(nList)/10)))==0
				fprintf('  Grid frame %d/%d (%.1f%%)\n', idxN, numel(nList), 100*idxN/numel(nList));
			end
		end
		close(vwGrid);
		close(figG);
		fprintf('SchizoVis: grid digit sweep video complete (%s).\n', gridVideoName);
	end
end

% -------- Figure export --------
if exportFigure
	try
		if ~isfolder(exportFolder); mkdir(exportFolder); end
		if VisPolar && ~isempty(fig) && isgraphics(fig)
			figPathFig = fullfile(exportFolder, sprintf('%s_polar.fig', FileName));
			figPathPng = fullfile(exportFolder, sprintf('%s_polar.png', FileName));
			savefig(fig, figPathFig);
			exportgraphics(fig, figPathPng, 'Resolution', 150);
			fprintf('SchizoVis: polar figure exported -> %s , %s\n', figPathFig, figPathPng);
		end
		if VisGrid
			% regenerate / show last grid number (static) for figure export
			schizoNumber = sVals(end);
			gridPathFig = fullfile(exportFolder, sprintf('%s_grid.fig', FileName));
			gridPathPng = fullfile(exportFolder, sprintf('%s_grid.png', FileName));
			GridDigitVis(schizoNumber, baseRadix, '', nList(end), 'FractionDigits', precisionOrder, 'Mode', mode, 'Colormap', parula(max(64,floor(baseRadix)*4)), 'Transparent', false, 'Show', true);
			try
				figGStatic = gcf; savefig(figGStatic, gridPathFig); exportgraphics(figGStatic, gridPathPng, 'Resolution', 150);
				fprintf('SchizoVis: grid figure exported -> %s , %s\n', gridPathFig, gridPathPng);
			catch ME2
				warning(ME2.identifier,'Grid figure export failed: %s', ME2.message);
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

% Optional Grid visualization (non-export) when VisGrid is enabled and exportNumber not used
if VisGrid
	try
		schizoNumber = sVals(end);
		% Show-only call (no file) - static grid view of last sqrt
		GridDigitVis(schizoNumber, baseRadix, '', nList(end), 'FractionDigits', precisionOrder, 'Mode', mode, 'Colormap', parula(max(64,floor(baseRadix)*4)), 'Transparent', true, 'Show', true);
		% Dock the created figure if any
		try
			figGrid = gcf;
			set(figGrid,'WindowStyle','docked');
		catch
			% ignore if figure cannot be docked
		end
	catch ME
		warning(ME.identifier, 'Could not create grid visualization: %s', ME.message);
	end
end