function [fVals, sVals] = SchizoGen(n, baseRadix)
% SchizoGen Build concatenated integer(s) 1..n (symbolic only)
% f(0)=0, f(k)=baseRadix*f(k-1)+k. Accepts scalar or vector n.
% Returns symbolic f(n) values and symbolic sqrt(f(n)); no vpa here.

arguments
    n double {mustBePositive}
    baseRadix double {mustBeGreaterThanOrEqual(baseRadix,2)} = 10 % allow non-integer (beta) bases >=2
end
if any(n~=floor(n)), error('n must be integer-valued'); end

% Treat baseRadix symbolically (may be non-integer). For non-integer bases this
% defines a beta-expansion style sequence f(k)=beta*f(k-1)+k.
vb = sym(baseRadix);
nVec = round(n(:));
maxN = max(nVec);
requested = unique(nVec);
markRequested = false(maxN,1); markRequested(requested) = true;

f = sym(0);
store = sym(zeros(numel(requested),1));
reqIdxMap = containers.Map(num2cell(requested), num2cell(1:numel(requested)));
for k=1:maxN
    f = vb*f + sym(k);
    if markRequested(k)
        store(reqIdxMap(k)) = f;
    end
end

% Map back to input order
fSym = sym(zeros(size(nVec)));
for i=1:numel(nVec)
    fSym(i) = store(reqIdxMap(nVec(i)));
end
sSym = sqrt(fSym); % purely symbolic sqrt

fVals = reshape(fSym, size(n));
sVals = reshape(sSym, size(n));
if isscalar(n)
    fVals = fVals(1); sVals = sVals(1);
end
end
