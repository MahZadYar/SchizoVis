function [fVals, sVals] = SchizoGen(n, baseRadix)
%SCHIZOGEN Build concatenated integer sequence f(n) and its square root.
%   [fVals, sVals] = SchizoGen(n, baseRadix) generates the "schizophrenic
%   number" sequence f(n) and its symbolic square root s(n) = sqrt(f(n)).
%   The sequence is defined by the recurrence relation:
%     f(0) = 0
%     f(k) = baseRadix * f(k-1) + k
%   This is equivalent to concatenating the digits of integers 1 through n
%   in a base specified by baseRadix.
%
%   The function uses symbolic computation to avoid premature rounding,
%   returning exact symbolic expressions for f(n) and its square root. It is
%   vectorized to handle multiple values of n efficiently.
%
%   Inputs:
%       n         - A scalar or vector of positive integers.
%       baseRadix - The numeral base (radix) to use for the sequence.
%                   Must be a numeric scalar greater than 1. Default is 10.
%
%   Outputs:
%       fVals     - A symbolic array of the same size as n, containing the
%                   values of f(n).
%       sVals     - A symbolic array of the same size as n, containing the
%                   symbolic square roots of f(n).
%
%   Example:
%       % Calculate f(3) and its square root in base 10
%       [f3, s3] = SchizoGen(3, 10);
%       % f3 will be the symbolic representation of 123
%       % s3 will be the symbolic representation of sqrt(123)
%
%       % Calculate for a vector of n values in base 2
%       n_vector = [1, 2, 3, 4];
%       [f_binary, s_binary] = SchizoGen(n_vector, 2);

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
