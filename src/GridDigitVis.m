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
%   'Mode'           : 'digits' (default) or 'residual' (continuous base expansion values)
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
opts.Mode = 'digits';
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
            case 'mode';           opts.Mode = char(lower(string(val)));
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

mode = lower(string(opts.Mode));
if ~(mode == "digits" || mode == "residual")
    warning('GridDigitViz:mode','Unknown Mode %s -> falling back to digits.', opts.Mode);
    mode = "digits";
end

% Generate value sequence according to mode
if mode == "digits"
    % Integer / fractional separation (beta base allowed) into discrete digits
    intPart = floor(symVal);
    fracPart = symVal - intPart;
    intDigits = [];
    if intPart == 0
        intDigits = 0;
    else
        tmp = intPart; loopCap = 5e5; it=0;
        while tmp > 0 && it < loopCap
            q = floor(tmp / baseRadix);
            r = tmp - q*baseRadix;
            d = double(floor(r));
            intDigits = [d intDigits]; %#ok<AGROW>
            tmp = q; it = it + 1;
        end
        if it >= loopCap
            warning('GridDigitViz:intLoopCap','Integer digit extraction reached loop cap; result may be incomplete.');
        end
    end
    fracDigits = [];
    if opts.FractionDigits > 0 && fracPart ~= 0
        decPrec = ceil(opts.FractionDigits * log(baseRadix)/log(10)) + 10;
        fnum = vpa(fracPart, decPrec);
        for k=1:opts.FractionDigits
            fnum = fnum * baseRadix;
            d = floor(fnum);
            fracDigits(end+1) = double(d); %#ok<AGROW>
            fnum = fnum - d;
        end
    end
    seqVals = [intDigits fracDigits];
else
    % Residual mode: capture residual values each step rather than discrete digits
    % We still perform an expansion but store the residual BEFORE extracting digit.
    intPart = floor(symVal);
    fracPart = symVal - intPart;
    seqVals = [];
    % integer residual sweep (divide backwards?) -> for visualization, build forward digits but store normalized residual r/baseRadix
    tmp = symVal;
    loopCap = 5e5; it=0; stackInt = [];
    while tmp > 0 && it < loopCap
        q = floor(tmp / baseRadix);
        r = tmp - q*baseRadix;
        % store residual normalized to base span
        stackInt = [double(r/baseRadix) stackInt]; %#ok<AGROW>
        tmp = q; it = it + 1;
    end
    if it >= loopCap
        warning('GridDigitViz:intLoopCap','Integer residual extraction loop cap reached.');
    end
    % fractional residual sweep
    if opts.FractionDigits > 0 && fracPart ~= 0
        decPrec = ceil(opts.FractionDigits * log(baseRadix)/log(10)) + 10;
        fnum = vpa(fracPart, decPrec);
        for k=1:opts.FractionDigits
            fnum = fnum * baseRadix;
            dWhole = floor(fnum);
            resNorm = double((fnum - dWhole)); % already fractional part after removing digit
            seqVals = [seqVals resNorm]; %#ok<AGROW>
            fnum = fnum - dWhole;
        end
    end
    seqVals = [stackInt seqVals];
    if isempty(seqVals); seqVals = 0; end
end

numDigits = numel(seqVals);
if numDigits == 0
    seqVals = 0; numDigits = 1;
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
    seqPadded = [seqVals padVals];
else
    seqPadded = seqVals;
end

img = reshape(seqPadded, cols, rows)'; % rows x cols

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
if mode == "digits"
    imgD(validMask) = max(0, min(digitCount-1, floor(imgD(validMask))));
    caxisRange = [0 digitCount-1];
else
    % residuals: values assumed in [0,1)
    % clamp to [0,1) then scale to actual base units (0 .. baseRadix)
    imgD(validMask) = max(0,min(1,imgD(validMask)));
    imgD(validMask) = imgD(validMask) * baseRadix;
    caxisRange = [0 baseRadix];
end

% Display image with NaN -> transparent; ensure first element maps to top-left
imAlpha = validMask;
imagesc(ax, imgD, 'AlphaData', imAlpha);
set(ax,'YDir','reverse'); % first row at top
axis(ax,'image'); axis(ax,'off');
if mode == "digits"
    colormap(ax, opts.Colormap(1:digitCount,:));
else
    % For residuals, ensure a sufficiently smooth map
    colormap(ax, interp1(linspace(0,1,size(opts.Colormap,1)), opts.Colormap, linspace(0,1,256)));
end
caxis(ax,caxisRange);
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
if mode == "digits"
    cb.Label.String = sprintf('Digits (beta %.4g, alphabet 0..%d)', baseRadix, digitCount-1);
else
    % Show residual scale in units of the base (0 .. baseRadix)
    cb.Label.String = sprintf('Residual (beta %.4g)', baseRadix);
    % choose 6 ticks across [0, baseRadix]
    tickVals = linspace(0, caxisRange(2), 6);
    cb.Ticks = tickVals;
    % format labels: integer-like for integer bases, else two decimals
    if baseRadix == floor(baseRadix)
        cb.TickLabels = compose('%d', round(tickVals));
    else
        cb.TickLabels = compose('%.2f', tickVals);
    end
end

if isempty(opts.Title)
    if mode == "digits"
        if opts.FractionDigits > 0
            ttl = sprintf('Schizo sqrt(f(%d)) beta %.4g digits (frac %d)', nVal, baseRadix, opts.FractionDigits);
        else
            ttl = sprintf('Schizo sqrt(f(%d)) beta %.4g digits', nVal, baseRadix);
        end
    else
        if opts.FractionDigits > 0
            ttl = sprintf('Schizo sqrt(f(%d)) beta %.4g residual (frac %d)', nVal, baseRadix, opts.FractionDigits);
        else
            ttl = sprintf('Schizo sqrt(f(%d)) beta %.4g residual', nVal, baseRadix);
        end
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
