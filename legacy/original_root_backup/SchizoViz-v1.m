%% Schizophrenic Number Path Visualizer
clear; clc;

%% User Configuration
nMin = 1;                  % Lowest n (we will use odd n only)
nMax = 301;                  % Highest n (we will use odd n only)
% Step size for n (if 2, will take every other n starting at nMin -> odd or even)
nStepSize = 2; % set to 1 for every n, 2 for alternating (odd/even based on nMin), etc.

% Ensure sane bounds
if nMin < 1
    warning('nMin must be >= 1. Clamping to 1.');
    nMin = 1;
end

if nMax < nMin
    error('nMax (%d) must be >= nMin (%d).', nMax, nMin);
end

% Ensure step is a positive integer
if nStepSize < 1 || mod(nStepSize,1)~=0
    error('nStepSize must be a positive integer.');
end

% Note: we keep nMin as the first sampled n. If you want to force the next
% odd/even value automatically, set nMin accordingly before running.
precisionOrder = 1500;      % Lowest exponent (e.g., 1200)
colormapName = 'parula';
reverseColormap = false;
export3d = true;    % enable 3D export
exportVideo = true;  % export sweep video of the final plot
export3Dplot = true; % export 3D plot of the final result
% Video filename includes n-range and precisionOrder for easy identification
videoFilename = sprintf('SchizoSweep_n%d-%d_p%d', nMin, nMax, precisionOrder);

videoFrameRate = 15;                % frames per second
videoView = 'top';                  % 'top' or '3d'
videoShowClock = true;              % draw digit clock per frame
videoMaxPointsPerFrame = inf;       % set to e.g. 20000 to subsample for speed
videoShowFullNumber = false;        % show full f(n) text (can be huge). If false, shortened.
videoLossless = true;               % use lossless 'Uncompressed AVI' when true (very large files)
videoResolution = [1080, 1080];      % [width height] in pixels; empty -> use current figure size

%% 1. Build f(n) values (only odd n) incrementally (symbolic, no precision yet)
fprintf('Building concatenated integer sequence up to n=%d...\n', nMax);
f_numbers = {};                  % cell array of sym values for odd n
% n values from nMin..nMax with user-specified step
nList = nMin:nStepSize:nMax;
f_n = sym(0);
idx = 0;
for k = 1:nMax
    f_n = f_n*10 + k;
    % Store f_n only for k values that match the sampling defined by nMin and nStepSize
    if k >= nMin && mod(k - nMin, nStepSize) == 0
        idx = idx + 1;
        f_numbers{idx} = f_n; %#ok<AGROW>
    end
end
fprintf('Stored %d n values (step=%d).\n', numel(f_numbers), nStepSize);

%% 2. Determine VPA precision from the largest f(n)
len_f_max = length(char(f_n));
digits_int_est = ceil(len_f_max/2);          % ~ integer digits of sqrt(f_nMax)
digits_frac_needed = abs(precisionOrder);
vpaDigits = digits_int_est + digits_frac_needed + 50;
digits(vpaDigits);
fprintf('VPA precision set to %d digits.\n', vpaDigits);

%% 3. Prepare colormap
% call colormap function by name
cmFunc = str2func(colormapName);
baseMap = cmFunc(256);
if reverseColormap
    baseMap = flipud(baseMap);
end

%% 4. Iterate odd n, extract per-exponent continuous slots, build points
% (per-layer data will be accumulated in matrices later; no need for early vectors)

p_min = -abs(precisionOrder);          % global lowest exponent
layerCount = numel(nList);

fprintf('Preparing uniform exponent range across layers...\n');
% Use largest n (last in list) to define global maximum exponent
s_max = vpa(sqrt(f_numbers{end}));
global_p_max = floor(double(log10(s_max)));
exponents = global_p_max:-1:p_min;  % uniform exponent vector for all layers
stepCountAll = numel(exponents);
fprintf('Uniform exponent range: [%d .. %d] (%d magnitudes)\n', global_p_max, p_min, stepCountAll);

% Precompute normalized z per layer in range [0,2]
if nMax == nMin
    zLayerVals = ones(1,layerCount) * 1; % single layer at center (1)
else
    zLayerVals = 2 * ( (nList - nMin) ./ (nMax - nMin) ); % maps nMin->0, nMax->2
end

% Radius normalization: map global_p_max -> r=0, exponent 0 -> r=1 (vector form)
if global_p_max > 0
    r_vals = 1 - (exponents ./ global_p_max);
else
    if p_min == 0
        r_vals = ones(size(exponents));
    else
        r_vals = (exponents - p_min) ./ (0 - p_min);
    end
end

% Allocate matrices: layers x magnitudes
digits_matrix = zeros(layerCount, stepCountAll);
xMat = zeros(layerCount, stepCountAll);
yMat = zeros(layerCount, stepCountAll);
zMat = zeros(layerCount, stepCountAll);
thetaMat = zeros(layerCount, stepCountAll);
colorMat = zeros(layerCount*stepCountAll, 3);

threshold_vpa = vpa(10)^p_min; % used only if we ever want to early-stop fractional tail (currently not used for full matrix)

for li = 1:layerCount
    n = nList(li);
    fprintf('Layer %d/%d (n=%d): computing sqrt & slots...\n', li, layerCount, n);
    s_vpa = vpa(sqrt(f_numbers{li}));
    residual = s_vpa;
    for jj = 1:stepCountAll
        e = exponents(jj);
        scale = vpa(10)^e;
        rawScaled = residual / scale;
        slot_vpa = mod(rawScaled, vpa(10));
        digits_matrix(li,jj) = double(slot_vpa);
        intPart = floor(rawScaled);
        residual = residual - intPart * scale; % will become zero; still continue to fill
    end
    theta_vals = (digits_matrix(li,:)/10) * 2*pi;
    xMat(li,:) = r_vals .* cos(theta_vals);
    yMat(li,:) = r_vals .* sin(theta_vals);
    thetaMat(li,:) = theta_vals;
    zMat(li,:) = zLayerVals(li); % broadcast
    % Color by layer height (z)
    zVal = zLayerVals(li);
    z01 = zVal / 2; % 0..1
    cmapIdx = 1 + z01 * (size(baseMap,1)-1);
    layerColor = interp1(1:size(baseMap,1), baseMap, cmapIdx);
    rowIdx = (li-1)*stepCountAll + (1:stepCountAll);
    colorMat(rowIdx,:) = repmat(layerColor, stepCountAll, 1);
end

% Flatten for scatter
% Flatten matrices in row-major (layer-major) order so colors align by layer
allX = reshape(xMat', 1, []);
allY = reshape(yMat', 1, []);
allZ = reshape(zMat', 1, []);
allC = colorMat; % rows already arranged per-layer
% allExp not used for plotting; omitted or compute if needed
allExp = reshape(repmat(exponents, layerCount, 1)', 1, []);
fprintf('Total points: %d (layers=%d * magnitudes=%d)\n', numel(allX), layerCount, stepCountAll);

%% 5. Plot stacked layers
figMain = figure('Name', sprintf('Stacked Prefix Map up to n=%d', nMax), 'Color',[0.08 0.08 0.1]);
ax = axes('Parent',figMain,'Color',[0.08 0.08 0.1]); hold(ax,'on'); set(figMain,'Renderer','opengl');
hScatter = scatter3(ax, allX, allY, allZ, 3, allC, 'filled', 'MarkerEdgeColor','none');

% Place the digit clock outside the data and draw repeating-digit guidelines
if isempty(allX)
    maxR = 1.0; else, maxR = max(sqrt(allX.^2 + allY.^2)); end
clockR = max(1.05, maxR * 1.12);
digit_angles = (0:9)/10 * 2*pi; tickInner = clockR * 0.97; tickOuter = clockR; labelR = clockR * 1.07;
for d=0:9
    ang = digit_angles(d+1);
    plot3(ax, [tickInner*cos(ang) tickOuter*cos(ang)], [tickInner*sin(ang) tickOuter*sin(ang)], [0 0], '-', 'Color',[0.5 0.5 0.5]);
    text(labelR*cos(ang), labelR*sin(ang), 0, sprintf('%d',d), 'Color','w','HorizontalAlignment','center','FontWeight','bold');
end
rep_angles = (1:9)/9 * 2*pi; guidR = clockR * 0.98;
for dd = 1:9
    ang = rep_angles(dd);
    plot3(ax, [0 guidR*cos(ang)], [0 guidR*sin(ang)], [0 0], '--', 'Color',[0.35 0.35 0.6], 'LineWidth', 0.8);
    text(guidR*1.06*cos(ang), guidR*1.06*sin(ang), 0, sprintf('%d.%s', dd, repmat(num2str(dd),1,3)), 'Color',[0.7 0.7 0.9], 'HorizontalAlignment','center');
end
thCirc = linspace(0,2*pi,360);
thCircR = 1;
plot3(ax, thCircR * cos(thCirc), thCircR * sin(thCirc), zeros(size(thCirc)), ':', 'Color',[0.6 0.6 0.6]);

axis(ax,'equal');
if ~isempty(allX)
    xmin = -clockR; xmax = clockR; xr = xmax - xmin; if xr == 0, xr = 1; end
    ymin = -clockR; ymax = clockR; yr = ymax - ymin; if yr == 0, yr = 1; end
    pad = 0.03; xlim(ax, [xmin - pad*xr, xmax + pad*xr]); ylim(ax, [ymin - pad*yr, ymax + pad*yr]);
else
    xlim(ax,[-1 1]); ylim(ax,[-1 1]);
end
if ~isempty(allZ)
    zmin = min(allZ); zmax = max(allZ);
    if zmin == zmax, m = 0.5; else, m = max((zmax - zmin) * 0.05, 0.05); end
    zlim(ax, [zmin - m, zmax + m]);
else
    zlim(ax,[0 2]);
end
title(ax, sprintf('Stacked Magnitude Sweep | odd n=%d..%d | precision [%d..%d]', nMin, nMax, 0, p_min), 'Color','w');
view(3);
grid(ax,'off'); ax.GridColor=[0.3 0.3 0.3]; ax.GridAlpha=0.25; ax.XColor='none'; ax.YColor='none'; ax.ZColor='none';
colormap(ax, baseMap);
hold(ax,'off');

% Optional: export the 3D figure as a .fig file using the same base name
if exist('export3Dplot','var') && export3Dplot
    try
        [p, nm, ext] = fileparts(videoFilename);
        if isempty(nm)
            nm = 'SchizoViz';
        end
        figPath = fullfile(pwd, [nm '.fig']);
        savefig(figMain, figPath);
        fprintf('Saved 3D figure to %s\n', figPath);
    catch ME
        warning(ME.identifier, '%s', ME.message);
    end
end

%% 6. Optional video export (reuse existing scatter3)
if exportVideo
    fprintf('Beginning video export to %s ...\n', videoFilename);
    % Choose writer based on lossless flag
    try
        if videoLossless
            % Uncompressed AVI is lossless but produces very large files
            warnMsg = sprintf('Video export in lossless mode (Uncompressed AVI) will produce very large files.\n');
            warning(warnMsg);
            % Ensure extension is .avi
            [~,~,ext] = fileparts(videoFilename);
            if ~strcmpi(ext, '.avi')
                videoFilename = [videoFilename '.avi'];
            end
            vw = VideoWriter(videoFilename, 'Uncompressed AVI');
        else
            % Prefer MPEG-4 (H.264) when available
            try
                vw = VideoWriter(videoFilename, 'MPEG-4');
            catch
                warning('MPEG-4 not available, falling back to Motion JPEG AVI');
                vw = VideoWriter(strrep(videoFilename,'.mp4','.avi'), 'Motion JPEG AVI');
            end
        end
    catch ME
        warning(ME.identifier, '%s', ME.message);
        error('Could not create VideoWriter.');
    end
    vw.FrameRate = videoFrameRate; open(vw);
    origX = get(hScatter,'XData'); origY = get(hScatter,'YData'); origZ = get(hScatter,'ZData');
    origTitle = get(get(ax,'Title'),'String'); [origAz, origEl] = view(ax);
    % Optionally resize figure to requested resolution for consistent frames
    origPos = get(figMain, 'Position'); % [left bottom width height]
    resized = false;
    if ~isempty(videoResolution) && numel(videoResolution) == 2
        try
            newPos = origPos; newPos(3:4) = videoResolution;
            set(figMain, 'Position', newPos);
            drawnow;
            resized = true;
        catch ME
            warning(ME.identifier, 'Could not resize figure to requested resolution: %s', ME.message);
        end
    end
    xFrame = nan(size(origX)); yFrame = nan(size(origY)); zFrame = nan(size(origZ));
    for li = 1:layerCount
        layerIdx = (li-1)*stepCountAll + 1 : li*stepCountAll;
        if isfinite(videoMaxPointsPerFrame) && numel(layerIdx) > videoMaxPointsPerFrame
            subIdxRel = round(linspace(1, numel(layerIdx), videoMaxPointsPerFrame));
            layerIdxUse = layerIdx(subIdxRel);
        else
            layerIdxUse = layerIdx;
        end
        xFrame(:)=nan; yFrame(:)=nan; zFrame(:)=nan;
        xFrame(layerIdxUse)=origX(layerIdxUse); yFrame(layerIdxUse)=origY(layerIdxUse); zFrame(layerIdxUse)=origZ(layerIdxUse);
        set(hScatter,'XData',xFrame,'YData',yFrame,'ZData',zFrame);
        fStrFull = char(f_numbers{li});
        if videoShowFullNumber
            fDisplay = fStrFull;
        else
            maxLen = 80;
            if length(fStrFull) > maxLen
                fDisplay = sprintf('%s...%s (len=%d)', fStrFull(1:30), fStrFull(end-30+1:end), length(fStrFull));
            else
                fDisplay = fStrFull;
            end
        end
        title(ax, sprintf('n = %d', nList(li)), 'Interpreter','none','Color','w');
        switch lower(videoView)
            case 'top'; view(ax,2); otherwise; view(ax,origAz,origEl);
        end
        drawnow limitrate; % Limit redraw rate
        writeVideo(vw, getframe(figMain));
        if mod(li,10)==0 || li==1 || li==layerCount
            fprintf('  Wrote frame %d/%d (%.1f%%)\n', li, layerCount, 100*li/layerCount);
        end
    end
    % restore original figure size & scatter
    set(hScatter,'XData',origX,'YData',origY,'ZData',origZ);
    if resized
        set(figMain, 'Position', origPos); drawnow;
    end
    title(ax, origTitle); view(ax, origAz, origEl); close(vw);
    fprintf('Video export complete.\n');
end

% Optional: launch volumetric viewers helper (disabled/commented by default)
% The helper creates a voxelized density, smooths it, renders an isosurface and
% (if available) a volshow volume. Uncomment the call below to run it.
% irrvis_viewers(allX, allY, allZ, allC, 'gridSize', 64, 'sigma', 1, 'isoLevel', 0.12);

return;