function outPath = GridDigitVis(symVal, baseRadix, outFile, nVal, varargin)
% GridDigitViz  Visualize a (possibly non-integer) symbolic number in base-R as a pixel image.
% Each digit (integer part + optional fractional part) becomes one pixel colored via a colormap.
% outPath = GridDigitViz(symVal, baseRadix, outFile, nVal, Name,Value...)
%
% Name-Value options:
%   'FractionDigits' : number of fractional base-R digits (default 0)
%   'Colormap'       : colormap name or Nx3 matrix (default parula(baseRadix))
%   'Transparent'    : true/false PNG transparency (default false)
%   'Title'          : custom title string (auto if empty)
%   'PadValue'       : value used for padding (NaN => transparent gaps) (default NaN)
%   'Square'         : true to force square layout (default true)
%   'Show'           : show figure (default false)
%   'MaxSide'        : target figure side in pixels (default 1000)
%
% NOTE: For very large integer parts this performs symbolic division; may be slow.

opts.FractionDigits = 0;
opts.Colormap = [];
opts.Transparent = false;
opts.Title = '';
opts.PadValue = NaN;
opts.Square = true;
opts.Show = false;
opts.MaxSide = 1000;
if ~isempty(varargin)
    for k=1:2:numel(varargin)
        key = lower(string(varargin{k})); val = varargin{k+1};
        switch key
            case 'fractiondigits'; opts.FractionDigits = max(0,floor(double(val)));
            case 'colormap';       opts.Colormap = val;
            case 'transparent';    opts.Transparent = logical(val);
            case 'title';          opts.Title = char(val);
            case 'padvalue';       opts.PadValue = val;
            case 'square';         opts.Square = logical(val);
            case 'show';           opts.Show = logical(val);
            case 'maxside';        opts.MaxSide = double(val);
            otherwise; warning('GridDigitViz:unknownOption','Unknown option %s', key);
        end
    end
end

% Resolve digit count (alphabet size) for integer / beta base
digitCount = floor(baseRadix);
if digitCount < 2
    error('GridDigitViz:digitCount','floor(baseRadix) must be >= 2');
end
% Resolve colormap
if isempty(opts.Colormap)
    opts.Colormap = parula(digitCount);
elseif ischar(opts.Colormap) || isstring(opts.Colormap)
    try
        cf = str2func(char(opts.Colormap));
    opts.Colormap = cf(max(digitCount,64));
    catch
    warning('GridDigitViz:colormap','Could not resolve colormap %s, using parula.', string(opts.Colormap));
        opts.Colormap = parula(digitCount);
    end
end
if size(opts.Colormap,2)~=3
    error('GridDigitViz:colormapShape','Colormap must be Nx3');
end
% Resample provided colormap to exactly baseRadix colors (spread across original map)
nC = size(opts.Colormap,1);
if nC == 1
    % single color -> replicate
    opts.Colormap = repmat(opts.Colormap, digitCount, 1);
else
    old = linspace(0,1,nC);
    new = linspace(0,1,digitCount);
    opts.Colormap = interp1(old, opts.Colormap, new);
end

% --- Integer digit extraction ---
% Integer / fractional separation (beta base allowed)
intPart = floor(symVal);
fracPart = symVal - intPart;

intDigits = [];
if intPart == 0
    intDigits = 0;
else
    tmp = intPart;
    % Use symbolic division; break on safety (loop cap)
    loopCap = 5e5; it=0;
    while tmp > 0 && it < loopCap
        q = floor(tmp / baseRadix);
        r = tmp - q*baseRadix; % remainder in [0, baseRadix)
        d = double(floor(r)); % digit in 0..floor(base)-1
        intDigits = [d intDigits]; %#ok<AGROW>
        tmp = q; it = it + 1;
    end
    if it >= loopCap
    warning('GridDigitViz:intLoopCap','Integer digit extraction reached loop cap; result may be incomplete.');
    end
end

% --- Fraction digits ---
fracDigits = [];
if opts.FractionDigits > 0 && fracPart ~= 0
    % numeric high precision
    decPrec = ceil(opts.FractionDigits * log(baseRadix)/log(10)) + 10;
    fnum = vpa(fracPart, decPrec);
    for k=1:opts.FractionDigits
    fnum = fnum * baseRadix;
    d = floor(fnum);
    fracDigits(end+1) = double(d); %#ok<AGROW>
    fnum = fnum - d;
    end
end

allDigits = [intDigits fracDigits];
numDigits = numel(allDigits);
if numDigits == 0
    allDigits = 0; numDigits=1;
end

% Layout
if opts.Square
    side = ceil(sqrt(numDigits));
    rows = side; cols = side;
else
    % rectangular closer to aspect 16:9
    aspect = 16/9;
    rows = ceil(sqrt(numDigits / aspect));
    cols = ceil(numDigits / rows);
end

pad = rows*cols - numDigits;
if pad > 0
    padVals = repmat(opts.PadValue,1,pad);
    allDigitsPadded = [allDigits padVals];
else
    allDigitsPadded = allDigits;
end

img = reshape(allDigitsPadded, cols, rows)'; % rows x cols

% Figure
vis = 'off'; if opts.Show, vis='on'; end
figColor = 'white'; if opts.Transparent, figColor='none'; end
fig = figure('Visible',vis,'Color',figColor,'Units','pixels','Position',[100 100 opts.MaxSide opts.MaxSide]);
ax = axes('Parent',fig); hold(ax,'on');
if opts.Transparent, set(ax,'Color','none'); end

% Prepare image data: ensure double, clamp valid digit range, and preserve NaN for padding
imgD = double(img);
% clamp digits to valid range; leave NaN as-is
validMask = ~isnan(imgD);
imgD(validMask) = max(0, min(digitCount-1, floor(imgD(validMask))));

% Display image with NaN -> transparent; ensure first element maps to top-left
imAlpha = validMask;
imagesc(ax, imgD, 'AlphaData', imAlpha);
set(ax,'YDir','reverse'); % first row at top
axis(ax,'image'); axis(ax,'off');
colormap(ax, opts.Colormap(1:digitCount,:));
caxis(ax,[0 digitCount-1]);
cb = colorbar(ax,'southoutside');
ticks = 0:digitCount-1;
cb.Ticks = ticks;
% tick labels (0-9 A B ...)
lbl = strings(size(ticks));
for i=1:numel(ticks)
    d = ticks(i);
    if d < 10
        lbl(i) = string(d);
    else
        lbl(i) = char('A' + (d-10));
    end
end
cb.TickLabels = cellstr(lbl);
cb.Label.String = sprintf('Digits (beta %.4g, alphabet 0..%d)', baseRadix, digitCount-1);

if isempty(opts.Title)
    if opts.FractionDigits > 0
    ttl = sprintf('Schizo number sqrt(f(%d)) beta %.4g (frac %d)', nVal, baseRadix, opts.FractionDigits);
    else
    ttl = sprintf('Schizo number sqrt(f(%d)) beta %.4g', nVal, baseRadix);
    end
else
    ttl = opts.Title;
end
title(ax, ttl,'Interpreter','none');

% Export
outPath = '';
if ~isempty(outFile)
    try
        [pdir,~,~] = fileparts(outFile);
        if ~isempty(pdir) && ~isfolder(pdir), mkdir(pdir); end
        if opts.Transparent
            try
                exportgraphics(fig,outFile,'BackgroundColor','none');
            catch
                exportgraphics(fig,outFile,'BackgroundColor','white');
            end
        else
            exportgraphics(fig,outFile,'BackgroundColor',figColor);
        end
        outPath = outFile;
    catch ME
    warning(ME.identifier,'GridDigitViz: export failed: %s', ME.message);
    end
end

% Backwards compatibility alias
function outPath = renderNumberImage(varargin)
outPath = GridDigitViz(varargin{:});
end

if ~opts.Show
    close(fig);
end
end
