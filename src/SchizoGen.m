function [fVals, sVals] = SchizoGen(n, baseRadix)
%SCHIZOGEN Generates terms of the "schizophrenic" number sequence.
%   [fVals, sVals] = SCHIZOGEN(n, baseRadix) computes the sequence defined
%   by f(0) = 0 and f(k) = baseRadix * f(k-1) + k for each integer k in the
%   input vector n. The function operates symbolically to avoid premature
%   rounding, making it suitable for high-precision analysis.
%
%   This process is equivalent to concatenating digits in the given base.
%   For example, in base 10, f(3) is 123. For non-integer bases, it
%   produces a generalized beta-expansion style sequence.
%
%   Syntax:
%       [fVals, sVals] = SchizoGen(n)
%       [fVals, sVals] = SchizoGen(n, baseRadix)
%
%   Input Arguments:
%       - n (double) - A scalar or vector of positive integers specifying
%         the sequence terms to compute.
%       - baseRadix (double) - The numeral base to use for the sequence
%         generation. Must be a value greater than 1. Defaults to 10.
%
%   Output Arguments:
%       - fVals (sym) - A symbolic array containing the computed f(k)
%         values, with the same dimensions as the input n.
%       - sVals (sym) - A symbolic array containing the square roots of the
%         f(k) values (sqrt(f(k))), with the same dimensions as n.
%
%   Examples:
%       % Generate the first 5 terms of the sequence in base 10
%       n_values = 1:5;
%       [f, s] = SchizoGen(n_values, 10);
%       disp('f(k) values:');
%       disp(f);
%       disp('sqrt(f(k)) values:');
%       disp(s);
%
%       % Generate a single term in a non-integer base
%       [f_beta, s_beta] = SchizoGen(3, 1.5);
%       disp(f_beta);
%
%   See also: ExpoExpand, PolarDigitVis, GridDigitVis

arguments
    n double {mustBePositive}
    baseRadix double {mustBeGreaterThan(baseRadix,1)} = 10 % allow non-integer (beta) bases >1
end
if any(n~=floor(n)), error('n must be integer-valued'); end

% Treat baseRadix symbolically (may be non-integer). For non-integer bases this
% defines a beta-expansion style sequence f(k)=beta*f(k-1)+k.
vb = sym(baseRadix);
nVec = round(n(:));
maxN = max(nVec);
requested = unique(nVec);
markRequested = false(maxN,1); markRequested(requested) = true;


% Vectorized Method:
exponents = (0:maxN-1)';   % exponents for vb^j, j=0..maxN-1
% build coefficients coeff(k,j)=max(0, k-j) so f(k)=sum_{j=0}^{k-1} coeff(k,j)*vb^j
% compute vb^j as symbolic column
% build coefficient matrix (numeric then convert to sym for reliable logical indexing)
coeff = (1:maxN)' - (0:maxN-1);  % maxN x maxN numeric by implicit expansion
coeff(coeff<0) = 0;
coeff = sym(coeff);             % convert to symbolic
cumTerms = coeff * vb.^exponents;  % maxN x 1 symbolic where cumTerms(k) = f(k)
fSym = cumTerms(nVec);          % pick values in input order


% % Itterative Method:
% f = sym(0);
% store = sym(zeros(numel(requested),1));
% reqIdxMap = containers.Map(num2cell(requested), num2cell(1:numel(requested)));
% for k=1:maxN
%     f = vb*f + sym(k);
%     if markRequested(k)
%         store(reqIdxMap(k)) = f;
%     end
% end

% % Map back to input order
% fSym = sym(zeros(size(nVec)));
% for i=1:numel(nVec)
%     fSym(i) = store(reqIdxMap(nVec(i)));
% end
% Use a symbolic half-power or the element-wise sqrt for symbolic arrays.
% The matrix-power operator ^ requires square matrices; here fSym is a vector
% (one f(k) per input n), so use element-wise operations. Using `sqrt` works
% for symbolic arrays and is clearer.
rootPower = sym(1)/2;       % symbolic 1/2 for explicitness
% element-wise symbolic power avoids using matrix ^ which requires square matrices
sSym = fSym .^ rootPower;   % per-element symbolic half-power

fVals = reshape(fSym, size(n));
sVals = reshape(sSym, size(n));
if isscalar(n)
    fVals = fVals(1); sVals = sVals(1);
end
