function D = expoExpandBatch(sVals, exponents, precisionOrder, baseRadix, mode)
%EXPOEXPANDBATCH Backward compatibility wrapper for ExpoExpand.
%   D = expoExpandBatch(sVals, exponents, precisionOrder, baseRadix, mode)
%   forwards all arguments to the `ExpoExpand` function.
%
%   This function is deprecated and is maintained for backward compatibility
%   with older scripts and tests. Use `ExpoExpand` in new code.
%
%   See also ExpoExpand.

if nargin < 5 || isempty(mode); mode = 'digits'; end
D = ExpoExpand(sVals, exponents, precisionOrder, baseRadix, mode);
end
