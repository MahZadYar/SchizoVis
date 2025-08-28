function [outPath, hStruct] = GridDigitVis(symVal, baseRadix, outFile, nVal, varargin)
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
%   'ColorScale'     : 'linear' (default) or 'exponential' - mapping for continuous residual colormap
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
opts.Mode = 'residual';
opts.ColorScale = 'exponential'; % 'linear' or 'exponential'
opts.Figure = [];
opts.FigureName = '';
opts.Tag = '';
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
            case 'colorscale';     opts.ColorScale = char(lower(string(val)));
            case 'figure';         opts.Figure = val;
            case 'figurename';     opts.FigureName = char(val);
            case 'tag';            opts.Tag = char(val);
            otherwise; warning('GridDigitViz:unknownOption','Unknown option %s', key);
        end
    end
end

% Resolve digit count (alphabet size) for integer / beta base
digitCount = ceil(baseRadix);
if digitCount < 1
    error('GridDigitViz:digitCount','floor(baseRadix) must be >= 1');
end
% Resolve colormap
if isempty(opts.Colormap)
    % Use full-resolution parula (256) to avoid low-base banding artifacts
    opts.Colormap = parula(256);
elseif ischar(opts.Colormap) || isstring(opts.Colormap)
    try
        cf = str2func(char(opts.Colormap));
        opts.Colormap = cf(256);
    catch
        warning('GridDigitViz:colormap','Could not resolve colormap %s, using parula(256).', string(opts.Colormap));
        opts.Colormap = parula(256);
    end
end
if size(opts.Colormap,2)~=3
    error('GridDigitViz:colormapShape','Colormap must be Nx3');
end
% Ensure we work with a high-resolution base map (>=256 rows recommended)
% Prepare mode and input vectorization
mode = string(opts.Mode);
multi = numel(symVal) > 1;
% ensure a column vector of values for iteration
symValVec = symVal(:);

% Nested helper: extract digit/residual sequence for a single value
    function seq = extractSeq(oneVal)
        % Support symbolic or numeric inputs
        isSym = isa(oneVal, 'sym');
        % Integer part
        intPart = floor(oneVal);
        % For 'digits' mode: produce integer digits (MSB->LSB) then fractional digits
        if mode == "digits"
            % Integer digits
            intDigits = [];
            if double(intPart) == 0
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
            % Fractional digits
            fracDigits = [];
            if opts.FractionDigits > 0
                fracPart = oneVal - floor(oneVal);
                if isSym
                    decPrec = ceil(opts.FractionDigits * log(baseRadix)/log(10)) + 10;
                    fnum = vpa(fracPart, decPrec);
                else
                    fnum = double(fracPart);
                end
                for kk=1:opts.FractionDigits
                    fnum = fnum * baseRadix;
                    dkk = floor(fnum);
                    fracDigits(end+1) = double(dkk); %#ok<AGROW>
                    fnum = fnum - dkk;
                end
            end
            seq = [intDigits fracDigits];
        else
            % 'residual' mode: produce integer residuals (normalized) then fractional residuals
            seq = [];
            % integer residuals (most significant first)
            tmp = oneVal; loopCap = 5e5; it=0; stackInt = [];
            while tmp > 0 && it < loopCap
                q = floor(tmp / baseRadix);
                r = tmp - q*baseRadix;
                stackInt = [double(r/baseRadix) stackInt]; %#ok<AGROW>
                tmp = q; it = it + 1;
            end
            if it >= loopCap
                warning('GridDigitViz:intLoopCap','Integer residual extraction loop cap reached.');
            end
            % fractional residuals
            if opts.FractionDigits > 0
                fracPart = oneVal - floor(oneVal);
                if isSym
                    decPrec = ceil(opts.FractionDigits * log(baseRadix)/log(10)) + 10;
                    fnum = vpa(fracPart, decPrec);
                else
                    fnum = double(fracPart);
                end
                for kk=1:opts.FractionDigits
                    fnum = fnum * baseRadix;
                    dWhole = floor(fnum);
                    resNorm = double((fnum - dWhole));
                    seq = [seq resNorm]; %#ok<AGROW>
                    fnum = fnum - dWhole;
                end
            end
            seq = [stackInt seq];
            if isempty(seq), seq = 0; end
        end
    end

if ~multi
    % Single value (retain original square packing behaviour)
    seqVals = extractSeq(symValVec);
    numDigits = numel(seqVals);
    if numDigits == 0, seqVals = 0; numDigits = 1; end
    if opts.Square
        side = ceil(sqrt(numDigits)); rows = side; cols = side;
    else
        aspect = 16/9; rows = ceil(sqrt(numDigits / aspect)); cols = ceil(numDigits / rows);
    end
    pad = rows*cols - numDigits;
    if pad > 0
        padVals = repmat(opts.PadValue,1,pad);
        seqPadded = [seqVals padVals];
    else
        seqPadded = seqVals;
    end
    img = reshape(seqPadded, cols, rows)';
else
    % Multi-n: produce one row per n (ignore square packing)
    nCount = numel(symValVec);
    seqCell = cell(nCount,1); maxLen = 0;
    for ii=1:nCount
        sRow = extractSeq(symValVec(ii));
        seqCell{ii} = sRow; maxLen = max(maxLen, numel(sRow));
    end
    if opts.Square
        % square doesn't apply in multi-row mode; warn once
        if nCount > 1
            warning('GridDigitViz:multiSquareIgnored','Square layout ignored for multi-row visualization.');
        end
    end
    % Align decimal point by padding integer part with leading zeros
    % Determine integer lengths per row (assume fractional tail length = opts.FractionDigits)
    intLens = zeros(nCount,1);
    for ii=1:nCount
        rlen = numel(seqCell{ii});
        fracLen = min(opts.FractionDigits, rlen);
        intLens(ii) = rlen - fracLen;
    end
    maxInt = max(intLens);
    maxLen = maxInt + max(0, opts.FractionDigits);
    img = nan(nCount, maxLen);
    for ii=1:nCount
        rowSeq = seqCell{ii};
        rlen = numel(rowSeq);
        fracLen = min(opts.FractionDigits, rlen);
        intLen = rlen - fracLen;
        % leading zeros to align integer part
        leadPad = maxInt - intLen;
        if leadPad > 0
            img(ii, 1:leadPad) = 0; % leading zeros
        end
        % place existing digits
        startIdx = leadPad + 1;
        img(ii, startIdx:startIdx + rlen -1) = rowSeq;
        % trailing pad with PadValue if shorter than full width
        if startIdx + rlen -1 < maxLen
            if ~isnan(opts.PadValue)
                img(ii, startIdx + rlen : maxLen) = opts.PadValue;
            else
                img(ii, startIdx + rlen : maxLen) = NaN;
            end
        end
    end
end

% Figure (allow reuse)
externalFig = ~isempty(opts.Figure) && isgraphics(opts.Figure,'figure');
vis = 'off'; if opts.Show, vis='on'; end
figColor = 'white'; if opts.Transparent, figColor='none'; end
if externalFig
    fig = opts.Figure;
    figure(fig);
    cla reset; % clear existing axes/graphics
    set(fig,'Visible',vis,'Color',figColor);
    if ~isempty(opts.FigureName)
        try
            set(fig,'Name',opts.FigureName);
        catch
        end
    end
    if ~isempty(opts.Tag)
        try
            set(fig,'Tag',opts.Tag);
        catch
        end
    end
else
    figArgs = {'Visible',vis,'Color',figColor,'Units','pixels','Position',[100 100 opts.MaxSide opts.MaxSide]};
    if ~isempty(opts.FigureName), figArgs = [figArgs {'Name',opts.FigureName}]; end %#ok<AGROW>
    if ~isempty(opts.Tag), figArgs = [figArgs {'Tag',opts.Tag}]; end %#ok<AGROW>
    fig = figure(figArgs{:});
end
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
imgH = imagesc(ax, imgD, 'AlphaData', imAlpha);
set(ax,'YDir','reverse'); % first row at top
% Fill available axes area (do not force equal axis) so image stretches to figure
axis(ax,'tight');
% leave axes position space for a right-side colorbar
try
    set(ax,'Units','normalized','Position',[0.06 0.07 0.86 0.88]);
catch
end

% If multi-row, show axis ticks: y -> n values, x -> exponent positions
if multi
    % show axes and label major ticks only
    axis(ax,'on');
    % y ticks: choose a few major ticks (up to 6) to label n values
    try
        nCountLocal = size(img,1);
        nMajor = min(6, nCountLocal);
        yTickPos = unique(round(linspace(1, nCountLocal, nMajor)));
        if isnumeric(nVal)
            yTickLabels = arrayfun(@(x) sprintf('%d',x), nVal(yTickPos), 'UniformOutput', false);
        else
            yTickLabels = arrayfun(@(x) char(x), nVal(yTickPos), 'UniformOutput', false);
        end
        set(ax,'YTick', yTickPos, 'YTickLabel', yTickLabels);
    catch
        set(ax,'YTick', unique(round(linspace(1,size(img,1),min(6,size(img,1))))));
    end
    % x ticks: show exponent numeric values (no leading 'p')
    try
        seqLens = cellfun(@numel, seqCell);
        intCounts = max(0, seqLens - opts.FractionDigits);
        maxInt = max(intCounts);
        leftExp = maxInt - 1;
        maxLenLocal = size(img,2);
        colExps = leftExp - (0:(maxLenLocal-1));
        nXT = min(8, maxLenLocal);
        xTickPos = unique(round(linspace(1, maxLenLocal, nXT)));
        xTickLabels = arrayfun(@(p) sprintf('%d',colExps(p)), xTickPos, 'UniformOutput', false);
        set(ax,'XTick', xTickPos, 'XTickLabel', xTickLabels);
        xlabel(ax,'exponent'); ylabel(ax,'n');
    catch
        set(ax,'XTick',[]);
    end
else
    axis(ax,'off');
end
% Choose colormap mapping depending on mode and ColorScale
if mode == "digits"
    % Discrete colormap for digits:
    % - if the provided colormap already has exactly digitCount entries
    %   (e.g. parula(digitCount)), use it as-is so the palette remains refined;
    % - otherwise, sample evenly across the provided colormap to produce
    %   digitCount colours (avoids truncating the palette to its first rows).
    if size(opts.Colormap,1) == digitCount
        map = opts.Colormap;
    else
        map = interp1(linspace(0,1,size(opts.Colormap,1)), opts.Colormap, linspace(0,1,digitCount));
    end
    colormap(ax, map);
else
    % For residuals, create a smooth colormap; support exponential scaling
    baseMap = opts.Colormap;
    xIn = linspace(0,1,size(baseMap,1));
    if strcmpi(opts.ColorScale,'exponential')
        xOut = (exp(linspace(0,1,256)) - 1) / (exp(1)-1); % exponential emphasis
    else
        xOut = linspace(0,1,256);
    end
    map = interp1(xIn, baseMap, xOut);
    colormap(ax, map);
end

caxis(ax,caxisRange);
cb = colorbar(ax,'eastoutside');
% Colorbar ticks: show only whole integers < baseRadix for digit labels
intMax = digitCount - 1;
% Always use whole-digit tick positions (0..digitCount-1) for the colorbar.
ticks = 0:intMax;
cb.Ticks = ticks;
if mode == "digits"
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
    cb.Label.String = sprintf('Digits (beta %.4g, alphabet 0..%d)', baseRadix, intMax);
else
    % For residuals, show the same integer tick positions (whole-digit markers)
    cb.TickLabels = compose('%d', ticks);
    cb.Label.String = sprintf('Residual (beta %.4g)', baseRadix);
end

% Build a plain title with root character and compact params
if isempty(opts.Title)
    if ~multi
        if mode == "digits"
            if opts.FractionDigits > 0
                ttl = sprintf('√ f(%d) | beta=%.4g | frac=%d | digits', nVal, baseRadix, opts.FractionDigits);
            else
                ttl = sprintf('√ f(%d) | beta=%.4g | digits', nVal, baseRadix);
            end
        else
            if opts.FractionDigits > 0
                ttl = sprintf('√ f(%d) | beta=%.4g | frac=%d | residual', nVal, baseRadix, opts.FractionDigits);
            else
                ttl = sprintf('√ f(%d) | beta=%.4g | residual', nVal, baseRadix);
            end
        end
    else
        nMinLocal = min(nVal); nMaxLocal = max(nVal);
        if mode == "digits"
            ttl = sprintf('√ f(n) | beta=%.4g | digits | n=[%d..%d]', baseRadix, nMinLocal, nMaxLocal);
        else
            ttl = sprintf('√ f(n) | beta=%.4g | residual | n=[%d..%d]', baseRadix, nMinLocal, nMaxLocal);
        end
    end
else
    ttl = opts.Title;
end
title(ax, ttl, 'Interpreter', 'none');

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
outPath = GridDigitVis(varargin{:});
end

if ~opts.Show && ~externalFig
    close(fig);
end

% Build handle struct if requested
if nargout > 1
    hStruct = struct('fig',fig,'ax',ax,'image',imgH,'colorbar',cb,'mode',char(mode), ...
        'baseRadix',baseRadix,'nVal',nVal,'multi',multi);
else
    hStruct = struct();
end
end
