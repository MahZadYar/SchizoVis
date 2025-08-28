function digitsOrResiduals = ExpoExpand(sVals, exponents, precisionOrder, baseRadix, mode)
% ExpoExpand  Vectorized expansion of multiple sqrt values across exponent range.
%   digitsOrResiduals = ExpoExpand(sVals, exponents, precisionOrder, baseRadix, mode)
%   sVals        : Lx1 (or 1xL) symbolic sqrt values (sqrt(f(n))) kept symbolic until here.
%   exponents    : 1xM descending exponent vector (same for all layers).
%   precisionOrder: nonnegative integer specifying negative exponent depth (already reflected in exponents).
%   baseRadix    : base > 1.
%   mode         : "digits" or "residual".
% Returns an LxM numeric matrix where each row corresponds to one sVal expansion.

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
