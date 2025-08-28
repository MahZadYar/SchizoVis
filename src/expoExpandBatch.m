function D = expoExpandBatch(sVals, exponents, precisionOrder, baseRadix, mode)
% expoExpandBatch  Backward compatibility wrapper for ExpoExpand (deprecated name).
%   D = expoExpandBatch(sVals, exponents, precisionOrder, baseRadix, mode)
%   forwards to ExpoExpand. Retained so legacy scripts/tests keep working.
%   Use ExpoExpand directly in new code.
if nargin < 5 || isempty(mode); mode = 'digits'; end
D = ExpoExpand(sVals, exponents, precisionOrder, baseRadix, mode);
end
