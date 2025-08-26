function h = PolarDigitVis(exponents, digitsMat, baseRadix, varargin)
% PolarDigitVis  Unified 2D/3D stacked radial digit plot (V2 style).
%   h = PolarDigitVis(E, D, base) where:
%     E : 1xM exponent vector (descending)
%     D : LxM digit/residual matrix OR 1xM vector (treated as single layer)
%     baseRadix : numeral base
%   Produces a clock-like radial layout per exponent column, stacking
%   layers along Z in [0, zMax].

if isvector(digitsMat)
    digitsMat = digitsMat(:)'; % force 1xM
end
exponents = exponents(:)';
[L,M] = size(digitsMat);
if numel(exponents)~=M
    error('Length mismatch: numel(exponents)=%d, size(digitsMat,2)=%d', numel(exponents), M);
end

% Parse options
p = inputParser;
p.addParameter('Figure', []);
p.addParameter('Colormap', parula(256));
p.addParameter('ReverseColormap', false, @(b)islogical(b)&&isscalar(b));
p.addParameter('Title', '', @(s)isstring(s)||ischar(s));
p.addParameter('LayerZRange', [], @(v)isnumeric(v)&&numel(v)==2&&v(1)<=v(2));
p.addParameter('DigitGuides', true, @(b)islogical(b)&&isscalar(b));
p.addParameter('RepeatGuides', true, @(b)islogical(b)&&isscalar(b));
p.addParameter('MarkerSize', 4, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('ColorMode', 'layer', @(s)any(strcmpi(s,{'layer','point'})));
p.addParameter('View', 'top', @(s) any(strcmpi(s,{'top','3d'})));
p.parse(varargin{:});
S = p.Results;

map = S.Colormap; if S.ReverseColormap, map = flipud(map); end

% Timing: total
t_total = tic;

% Figure / axes
if isempty(S.Figure)
    fig = figure('Color',[0.08 0.08 0.1], 'Name','PolarDigitVis Stack');
    % Prefer OpenGL for better scatter3 performance and proper rendering
    try set(fig,'Renderer','opengl'); catch, end
else
    fig = S.Figure; figure(fig);
end
ax = axes('Parent', fig, 'Color',[0.08 0.08 0.1]); hold(ax,'on');

% Radial mapping: max exponent -> radius 0; exponent 0 -> 1
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

% Angles (vectorized). In residual mode values still in [0, baseRadix)
t_coords = tic;
theta = (digitsMat / baseRadix) * 2*pi; % LxM
% Broadcast radii to layers
rMat = repmat(r_vals, L,1);
xMat = cos(theta) .* rMat;
yMat = sin(theta) .* rMat;
coordsTime = toc(t_coords);
fprintf('PolarDigitVis: coords build %.3fs (%d points)\n', coordsTime, numel(xMat));

% Clock radius: use maximum mapped radius so labels/guides scale with data
clockR = max(r_vals);

% Compute zRange defaulted to [0 clockR] if user left it empty
if isempty(S.LayerZRange)
    zRange = [0, clockR];
else
    zRange = S.LayerZRange;
end
zVals = linspace(zRange(1), zRange(2), max(L,2));
if L==1, zVals = mean(zRange); end
zMat = zVals(:) * ones(1,M);

% Colors: build LxM color matrix then flatten consistently
% Colors computation
t_color = tic;
switch lower(S.ColorMode)
    case 'layer'
        idxLayer = linspace(1, size(map,1), L); % one color per layer
        layerCols = interp1(1:size(map,1), map, idxLayer); % Lx3
        colorMat = zeros(L,M,3);
        for c = 1:3
            colorMat(:,:,c) = repmat(layerCols(:,c), 1, M);
        end
    case 'point'
        % gradient over radius per column
        rNorm = (r_vals - min(r_vals))/(max(r_vals)-min(r_vals)+eps);
        idxPoint = 1 + rNorm*(size(map,1)-1);
        pointCols = interp1(1:size(map,1), map, idxPoint); % Mx3
        colorMat = zeros(L,M,3);
        for c = 1:3
            colorMat(:,:,c) = repmat(pointCols(:,c)', L, 1);
        end
    otherwise
        colorMat = repmat(reshape(map(1,:),1,1,3), L, M);
end
colorTime = toc(t_color);
fprintf('PolarDigitVis: color build %.3fs\n', colorTime);

% Flatten in same order as xMat'(:)
t_plot = tic;
allX = reshape(xMat', [], 1);
allY = reshape(yMat', [], 1);
allZ = reshape(zMat', [], 1);
% Build flattened color matrix aligned with xMat'(:) ordering
numPts = numel(allX);
allC = zeros(numPts, 3);
for c = 1:3
    allC(:,c) = reshape(colorMat(:,:,c)', [], 1);
end

% scatter3 plot
h = scatter3(ax, allX, allY, allZ, S.MarkerSize, allC, 'filled', 'MarkerEdgeColor','none');
% Disable hit testing to reduce UI overhead
try set(h,'HitTest','off'); catch, end
plotTime = toc(t_plot);
fprintf('PolarDigitVis: flatten+plot %.3fs\n', plotTime);

% Guides: outer dotted circle at unit radius (decimal point indicator)
thCirc = linspace(0,2*pi,360);
plot3(ax, cos(thCirc), sin(thCirc), ones(size(thCirc))*zRange(1), ':', 'Color',[0.6 0.6 0.6]);
% Clock outer circle at clockR
plot3(ax, clockR*cos(thCirc), clockR*sin(thCirc), ones(size(thCirc))*zRange(1), '-', 'Color',[0.6 0.6 0.6]);

% Digit spokes & labels
if S.DigitGuides
    whole_angles = (0:baseRadix-1)/baseRadix * 2*pi;
    for d = 0:baseRadix-1
        ang = whole_angles(d+1);
        plot3(ax, [0 clockR*cos(ang)], [0 clockR*sin(ang)], [zRange(1) zRange(1)], '-', 'Color',[0.5 0.5 0.5], 'LineWidth', 0.6);
        if d < 10, ds = num2str(d); else, ds = char('A'+(d-10)); end
        text(clockR*1.05*cos(ang), clockR*1.05*sin(ang), zRange(1), ds, 'Color',[0.9 0.9 0.9], 'HorizontalAlignment','center');
    end
end

% Simplified, vectorized guides to reduce plot load
% - sample radii if there are many layers
% - draw all major ticks in a single plot call (with NaN separators)
% - limit repeat guides to a small number and vectorize their ticks too

% sample radii to at most 60 ticks total per spoke
maxTicksPerSpoke = 60;
nR = numel(r_vals);
if nR > maxTicksPerSpoke
    idxR = unique(round(linspace(1, nR, maxTicksPerSpoke)));
    r_sample = r_vals(idxR);
else
    r_sample = r_vals;
end

tickLen = clockR * 0.02;
% build arrays with NaNs to plot all major ticks at once
Xticks = [];
Yticks = [];
Zticks = [];
for d = 0:baseRadix-1
    ang = (d/baseRadix) * 2*pi;
    cang = cos(ang); sang = sin(ang);
    for r = r_sample
        r_inner = max(0, r - tickLen);
        Xticks = [Xticks; r_inner*cang; r*cang; NaN];
        Yticks = [Yticks; r_inner*sang; r*sang; NaN];
        Zticks = [Zticks; zRange(1); zRange(1); NaN];
    end
end
% single call to plot all major ticks
if ~isempty(Xticks)
    plot3(ax, Xticks(:), Yticks(:), Zticks(:), '-', 'Color',[0.55 0.55 0.55], 'LineWidth', 0.7);
end

% Repetition guides: limit number and vectorize
if S.RepeatGuides && baseRadix>2
    repBase = baseRadix-1;
    maxRep = inf; % clamp to reasonable number of repetition guides
    nRep = min(repBase, maxRep);
    repIdx = unique(round(linspace(1, repBase, nRep)));
    rep_angles = (repIdx ./ repBase) * 2*pi;
    guidR = clockR * 0.98;
    % plot all repeat spokes in one call
    XR = [];
    YR = [];
    ZR = [];
    for ang = rep_angles
        XR = [XR; 0; guidR*cos(ang); NaN];
        YR = [YR; 0; guidR*sin(ang); NaN];
        ZR = [ZR; zRange(1); zRange(1); NaN];
    end
    if ~isempty(XR)
        plot3(ax, XR(:), YR(:), ZR(:), '--', 'Color',[0.35 0.35 0.6], 'LineWidth', 0.6);
    end

    % dashed ticks along these repeated spokes (also vectorized)
    tickLen2 = clockR * 0.015;
    XRt = [];
    YRt = [];
    ZRt = [];
    for ang = rep_angles
        cang = cos(ang); sang = sin(ang);
        for r = r_sample
            r_inner = max(0, r - tickLen2);
            XRt = [XRt; r_inner*cang; r*cang; NaN];
            YRt = [YRt; r_inner*sang; r*sang; NaN];
            ZRt = [ZRt; zRange(1); zRange(1); NaN];
        end
    end
    if ~isempty(XRt)
        plot3(ax, XRt(:), YRt(:), ZRt(:), '--', 'Color',[0.35 0.45 0.7], 'LineWidth', 0.6);
    end

    % sparse labels for repetition guides (only those plotted)
    for k = 1:numel(repIdx)
        dd = repIdx(k);
        ang = rep_angles(k);
        if dd < 10, symDig = num2str(dd); else, symDig = char('A'+(dd-10)); end
        repChar = char(8230); % ellipsis
        labelStr = sprintf('%s.%s', symDig, [symDig repChar]);
        text(guidR*1.03*cos(ang), guidR*1.03*sin(ang), zRange(1), labelStr, ...
            'Color',[0.7 0.7 0.9], 'HorizontalAlignment','center', 'FontSize', 9);
    end
end


axis(ax,'equal'); axis(ax,'off'); colormap(ax, map);
if strlength(S.Title)>0, title(ax, S.Title, 'Color','w'); end
% Ensure plot renders before reporting timings
drawnow;
totalTime = toc(t_total);
fprintf('PolarDigitVis: total %.3fs\n', totalTime);
% Apply requested view (default 'top')
switch lower(string(S.View))
    case 'top'
        view(ax,2);
    otherwise
        view(ax,3);
end
hold(ax,'off');
end
