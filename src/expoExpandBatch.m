function D = expoExpandBatch(sVals, exponents, precisionOrder, baseRadix, mode)
%EXPOEXPANDBATCH Backward compatibility wrapper for ExpoExpand.
%   D = EXPOEXPANDBATCH(...) is a deprecated function that forwards all its
%   arguments to ExpoExpand. It is retained for backward compatibility
%   with older scripts and tests.
%
%   Please use ExpoExpand directly in all new code.
%
%   See also: ExpoExpand
if nargin < 5 || isempty(mode); mode = 'digits'; end
D = ExpoExpand(sVals, exponents, precisionOrder, baseRadix, mode);
end
