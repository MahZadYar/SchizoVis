classdef testExpoExpand < matlab.unittest.TestCase
    methods(Test)
        function basicDigits(testCase)
            sVal = sym('sqrt(2)');
            exponents = [1 0 -1 -2];
            precisionOrder = 3;
            base = 10;
            mode = 'digits';
            D = ExpoExpand(sVal, exponents, precisionOrder, base, mode);
            testCase.verifySize(D,[1 numel(exponents)]);
            testCase.verifyGreaterThanOrEqual(D,0);
        end
    end
end
