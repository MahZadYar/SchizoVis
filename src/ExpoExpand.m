function digitsOrResiduals = ExpoExpand(sVals, exponents, precisionOrder, baseRadix, mode)
%EXPOEXPAND Vectorized base-b expansion of multiple symbolic values.
%   D = EXPOEXPAND(S, E, P, B, M) computes the base-b expansion for a
%   vector of symbolic values S across a shared grid of exponents E. This
%   function is optimized for expanding multiple numbers (e.g., sqrt(f(n))
%   from SchizoGen) simultaneously, which is more efficient than processing
%   them individually.
%
%   The function can operate in two modes: 'digits' mode, which extracts
%   the integer digits of the expansion, or 'residual' mode, which returns
%   the continuous residual values at each exponent step.
%
%   Syntax:
%       digitsOrResiduals = ExpoExpand(sVals, exponents, precisionOrder, baseRadix, mode)
%
%   Input Arguments:
%       - sVals (sym) - A vector (Lx1 or 1xL) of symbolic values to be
%         expanded. These are typically the sqrt(f(n)) values from SchizoGen.
%       - exponents (double) - A 1xM vector of descending exponents that
%         defines the sampling positions for the expansion.
%       - precisionOrder (double) - A non-negative integer that specifies
%         the depth of fractional exponents. This value is used for context
%         but is assumed to be already incorporated into the 'exponents' vector.
%       - baseRadix (double) - The numeral base (b > 1) for the expansion.
%         Can be a non-integer (for beta-expansions).
%       - mode (string) - The expansion mode. Must be either "digits" to
%         extract integer digits (0..floor(b)-1) or "residual" to get
%         continuous values in [0, b).
%
%   Output Arguments:
%       - digitsOrResiduals (double) - An LxM matrix where each row contains
%         the expansion of the corresponding value from sVals. The content
%         of the matrix depends on the selected 'mode'.
%
%   Notes:
%       The function includes an optional internal validation mechanism to
%       check the precision of the expansion. To enable it, set the
%       'enableValidation' flag to true within the function's body.
%
%   Examples:
%       % Generate symbolic square roots
%       [~, sVals] = SchizoGen(1:3, 10);
%       % Define an exponent range
%       exponents = 2:-1:-2;
%       precisionOrder = 2;
%       % Get integer digits
%       digitsMat = ExpoExpand(sVals, exponents, precisionOrder, 10, "digits");
%       disp(digitsMat);
%
%   See also: SchizoGen, PolarDigitVis, GridDigitVis, expoExpandBatch

arguments
    sVals {mustBeNonempty}
    exponents double {mustBeNonempty}
    precisionOrder double {mustBeNonnegative, mustBeInteger} %#ok<INUSA>
    baseRadix double {mustBeGreaterThan(baseRadix,1)} % allow non-integer base (beta)
    mode {mustBeTextScalar}
end

mode = string(lower(mode));
if ~(mode=="digits" || mode=="residual")
    error('mode must be "digits" or "residual"');
end

% Normalize shapes
exponents = exponents(:)';
L = numel(sVals);
M = numel(exponents);

% Use the global symbolic precision (set by digits(n) in the driver).
currDigits = digits; %#ok<DIGIT>
fprintf('ExpoExpand: using global symbolic precision digits=%d (driver-controlled, beta-base)\n', currDigits);

% Optional validation (lightweight). Set enableValidation=true to compare a single row with higher precision.
enableValidation = false; %#ok<NASGU> (user can toggle manually)
validationRow = 1;          % which row to test
validationExtra = 80;       % extra decimal digits for validation run
validationTail = 25;        % how many deepest (most negative exponent) slots to compare


% High-precision numeric conversion for all sqrt values at once using the
% global symbolic precision (no explicit second arg to vpa so current digits apply).
sNums = vpa(sVals(:)); % Lx1 symbolic high-precision using current digits
basePowMax = vpa(baseRadix)^exponents(1);
residuals = sNums / basePowMax; % Lx1 numeric symbolic

% Preallocate result
digitsOrResiduals = zeros(L, M);

for j = 1:M
    intParts = floor(residuals); % Lx1 symbolic ints
    if mode == "digits"
        % For non-integer base (beta-expansion) digit alphabet is 0..floor(beta)-1
        digitsOrResiduals(:,j) = double(intParts);
    else
        digitsOrResiduals(:,j) = double(residuals);
    end
    residuals = (residuals - intParts) * baseRadix;
end

% ---------------- Optional precision validation ----------------
if exist('enableValidation','var') && enableValidation && validationRow <= L
    try
        fprintf('ExpoExpand: validating precision on row %d ...\n', validationRow);
        sTest = sVals(validationRow);
    % Use a higher-precision copy for validation: current digits + extra
    currDigits = digits; %#ok<DIGIT>
    sNumHi = vpa(sTest, currDigits + validationExtra);
        resHi = sNumHi / basePowMax;
        valsHi = zeros(1,M);
        for j = 1:M
            intHi = floor(resHi);
            if mode == "digits"
                valsHi(j) = double(intHi);
            else
                valsHi(j) = double(resHi);
            end
            resHi = (resHi - intHi) * baseRadix;
        end
        tailIdx = max(1, M - validationTail + 1):M;
        baseVals = digitsOrResiduals(validationRow, tailIdx);
        hiVals = valsHi(tailIdx);
        diffVals = abs(baseVals - hiVals);
        maxDiff = max(diffVals);
        if mode == "digits"
            if any(diffVals > 0)
                warning('ExpoExpand: digit mismatch detected in validation tail (max diff %g). Consider increasing precision.', maxDiff);
            else
                fprintf('ExpoExpand: validation OK (no digit differences in tail of %d positions).\n', numel(tailIdx));
            end
        else
            tol = 1e-12; % residual numeric comparison tolerance
            if maxDiff > tol
                warning('ExpoExpand: residual mismatch (max diff %.3g > tol %.1g). Increase precision.', maxDiff, tol);
            else
                fprintf('ExpoExpand: residual validation OK (max diff %.3g <= tol %.1g).\n', maxDiff, tol);
            end
        end
    catch ME
        warning(ME.identifier,'Validation phase skipped: %s', ME.message);
    end
end

end
