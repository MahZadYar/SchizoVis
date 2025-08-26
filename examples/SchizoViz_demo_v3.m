%% SchizoViz v3 Demonstration Script
% Maintainer: Maziar Moussavi

clear; clc; close all;
nMin = 1; nMax = 101; nStepSize = 2; baseRadix = 12; precisionOrder = 1500; mode = "digits";
exportFigure=false; exportVideo=false; exportNumber=true; exportGridVideo=false; %#ok<NASGU>
videoLayers=inf; videoFrameRate=15; videoView='top'; videoLossless=true; videoResolution=[1080 1080]; videoLastN=inf; videoDynamicColors=true; %#ok<NASGU>
exportFolder = fullfile(pwd,'outputs'); if ~isfolder(exportFolder); mkdir(exportFolder); end
FileName = sprintf('SchizoViz_n%d-%d_b%d_p%d', nMin, nMax, baseRadix, precisionOrder);
nList = nMin:nStepSize:nMax; [fVals, sVals] = SchizoGen(nList, baseRadix); %#ok<NASGU>
p_min = -abs(precisionOrder); s_max_sym = sVals(end); global_p_max = ceil(double(log(vpa(s_max_sym))/log(baseRadix))); exponents = global_p_max:-1:p_min;
digitsMat = ExpoExpand(sVals, exponents, precisionOrder, baseRadix, mode);
PolarDigitVis(exponents, digitsMat, baseRadix, 'Title', sprintf('n=%d..%d base=%d', nMin, nMax, baseRadix));
if exportNumber
    fractionDigits = precisionOrder; schizoNumber = sVals(end); outFile = fullfile(exportFolder, sprintf('%s_schizo_n%d_base%d.png', FileName, nList(end), baseRadix)); cmapForNum = parula(baseRadix); %#ok<NASGU>
    GridDigitViz(schizoNumber, baseRadix, outFile, nList(end), 'FractionDigits', fractionDigits, 'Colormap', cmapForNum, 'Transparent', true, 'Show', false);
    fprintf('Exported digit grid -> %s\n', outFile);
end
fprintf('Demo complete.\n');
