function [fVals, sVals] = SchizoGen(n, baseRadix)
% SchizoGen Build concatenated integer(s) 1..n (symbolic only)
% f(0)=0, f(k)=baseRadix*f(k-1)+k. Accepts scalar or vector n.
% Returns symbolic f(n) values and symbolic sqrt(f(n)); no vpa here.

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


% % Classic Method:
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

sSym = sqrt(fSym); % purely symbolic sqrt

fVals = reshape(fSym, size(n));
sVals = reshape(sSym, size(n));
if isscalar(n)
    fVals = fVals(1); sVals = sVals(1);
end
end
