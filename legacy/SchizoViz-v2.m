%% Irrational (sqrt sequence) radial layer plot (Legacy v2)
% Original legacy script retained for historical reference (unmodified logic).
% For modern usage see ../examples/SchizoViz_demo_v3.m

%% ORIGINAL CONTENT BELOW
%% Irrational (sqrt sequence) radial layer plot
clear; clc;

% ---------------- Parameters ----------------
nMin = 1; nMax = 301;      % inclusive range of n (sampled with nStepSize)
nStepSize = 2;               % stride over n
precisionOrder = 2000;       % negative lowest exponent (precision depth)
baseRadix = 12;              % numeral base
colormapName = 'parula'; reverseColormap = false;
exportVideo = false; export3Dplot = false;
videoFrameRate = 15; videoView = 'top';
videoMaxPointsPerFrame = inf; videoShowFullNumber = false;
videoLossless = true; videoResolution = [1080 1080];
useGPU = true; continuousAngleRefinement = false; %#ok<NASGU> (kept for future)
videoFilename = sprintf('SchizoSweep_base%d_n%d-%d_p%d', baseRadix, nMin, nMax, precisionOrder);

% --------------- Validation -----------------
if baseRadix < 2 || baseRadix ~= floor(baseRadix), error('baseRadix integer >=2'); end
if nMin < 1, nMin = 1; end
if nMax < nMin, error('nMax < nMin'); end
if nStepSize < 1 || mod(nStepSize,1), error('nStepSize positive int'); end

%% 1. Build concatenated integer sequence f(n)
t_total = tic; t_build = tic;
nList = nMin:nStepSize:nMax; layerCount = numel(nList);
f_numbers = cell(1, layerCount);
f_n = sym(0); idx = 0;
for k = 1:nMax
    % Decimal-style concatenation: multiply by 10^digits(k) then add k
    f_n = f_n * baseRadix + k;     % append
    if k >= nMin && mod(k - nMin, nStepSize)==0
        idx = idx + 1;
        f_numbers{idx} = f_n; %#ok<AGROW>
    end
end
fprintf('Built %d sampled n values (stride=%d) in %.2fs\n', layerCount, nStepSize, toc(t_build));

%% 2. Precision (rough: integer digits of sqrt ~ half the digits of f_n)
t_vpa = tic;
len_f = length(char(f_n));
vpaDigits = ceil(len_f/2) + abs(precisionOrder) + 2;
digits(vpaDigits);
fprintf('VPA digits = %d (len(f_n)=%d) in %.2fs\n', vpaDigits, len_f, toc(t_vpa));

%% 3. Colormap
baseMap = feval(colormapName, 256); if reverseColormap, baseMap = flipud(baseMap); end

%% 4. Exponent grid (shared across layers)
p_min = -abs(precisionOrder);
s_max = vpa(sqrt(f_numbers{end}));
global_p_max = ceil(double(log(s_max)/log(baseRadix)));
exponents = global_p_max:-1:p_min; stepCountAll = numel(exponents);
fprintf('Exponent range [%d..%d] (%d steps)\n', global_p_max, p_min, stepCountAll);

% Layer z in [0,2]
if nMax == nMin, zLayerVals = ones(1,layerCount); else, zLayerVals = 2*(nList - nMin)/(nMax - nMin); end
digits_matrix = zeros(layerCount, stepCountAll);
if global_p_max > 0
    r_vals = 1 - exponents/global_p_max;
else
    r_vals = (exponents - p_min)/(0 - p_min + (p_min==0));
end

%% 5. Vectorized digit extraction
vb = vpa(baseRadix); basePowMax = vb^global_p_max;
t_extract = tic;
s_arr = vpa(sqrt([f_numbers{:}]));            % 1 x L
residuals = s_arr / basePowMax;
tenPct = max(1,floor(stepCountAll/10));
for j = 1:stepCountAll
    slot = mod(residuals, vb);
    digits_matrix(:,j) = double(slot(:));
    residuals = slot * baseRadix;
    if j==1 || j==stepCountAll || mod(j,tenPct)==0, fprintf('  exp %d/%d\n', j, stepCountAll); end
end
fprintf('Extraction %.2fs\n', toc(t_extract));

%% 6. Coordinate mapping
thetaMat = digits_matrix/baseRadix * (2*pi);
xMat = cos(thetaMat) .* r_vals; yMat = sin(thetaMat) .* r_vals; %#ok<NASGU>
zMat = zLayerVals(:) .* ones(1, stepCountAll);
% Colors per layer (single color each layer)
idxFloat = 1 + (zLayerVals/2) * (size(baseMap,1)-1);
colorsPerLayer = interp1(1:size(baseMap,1), baseMap, idxFloat);
colorMat = kron(colorsPerLayer, ones(stepCountAll,1));

%% 7. Flatten arrays
t_flat = tic;
allX = reshape(xMat',1,[]); allY = reshape(yMat',1,[]); allZ = reshape(zMat',1,[]);
allC = colorMat; allExp = reshape(repmat(exponents, layerCount,1)',1,[]); %#ok<NASGU>
fprintf('Points: %d (%dx%d) flatten %.2fs\n', numel(allX), layerCount, stepCountAll, toc(t_flat));
fprintf('Total so far %.2fs\n', toc(t_total));

%% 8. Plot
figMain = figure('Name', sprintf('Stacked Magnitude Sweep n=%d..%d', nMin, nMax), 'Color',[0.08 0.08 0.1]);
ax = axes('Parent',figMain,'Color',[0.08 0.08 0.1]); hold(ax,'on');
clockR = max(r_vals);
digit_angles = (0:baseRadix-1)/baseRadix * 2*pi; %#ok<NASGU>
hScatter = scatter3(ax, allX, allY, allZ, 3, allC, 'filled', 'MarkerEdgeColor','none'); %#ok<NASGU>

% Repetition guides (skip 0)
repBase = baseRadix-1;
rep_angles = (1:repBase)/repBase * 2*pi; guidR = clockR * 0.98;
for dd = 1:repBase
    ang = rep_angles(dd);
    plot3(ax, [0 guidR*cos(ang)], [0 guidR*sin(ang)], [0 0], '--', 'Color',[0.35 0.35 0.6], 'LineWidth', 0.8);
    % label symbol
    if dd < 10
        symDig = num2str(dd);
    else
        symDig = char('A' + (dd - 10));
    end
    repChar = char(8230);
    labelStr = sprintf('%s.%s', symDig, [symDig repChar]);
    text(guidR*1.06*cos(ang), guidR*1.06*sin(ang), 0, labelStr, 'Color',[0.7 0.7 0.9], 'HorizontalAlignment','center');
end

% Whole-digit radial spokes
whole_angles = (0:baseRadix-1)/baseRadix * 2*pi; wholeR = clockR * 1.02;
for dd = 0:(baseRadix-1)
    ang = whole_angles(dd+1);
    plot3(ax, [0 wholeR*cos(ang)], [0 wholeR*sin(ang)], [0 0], '-', 'Color',[0.5 0.5 0.5], 'LineWidth', 0.6);
    % label symbol
    if dd < 10
        symDig = num2str(dd);
    else
        symDig = char('A' + (dd - 10));
    end
    text(wholeR*1.08*cos(ang), wholeR*1.08*sin(ang), 0, symDig, 'Color',[0.9 0.9 0.9], 'HorizontalAlignment','center');
end
% Outer circle
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
title(ax, sprintf('sqrt(concat 1..n) base %d  n=%d..%d  exp[%d..%d]', baseRadix, nMin, nMax, global_p_max, p_min), 'Color','w');
view(3);
grid(ax,'off'); ax.GridColor=[0.3 0.3 0.3]; ax.GridAlpha=0.25; ax.XColor='none'; ax.YColor='none'; ax.ZColor='none';
colormap(ax, baseMap);
hold(ax,'off');

return;
