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
        function residualModeBasic(testCase)
            sVal = sym('sqrt(2)');
            exponents = 0:-1:-4;
            precisionOrder = 5;
            base = 10;
            mode = 'residual';
            R = ExpoExpand(sVal, exponents, precisionOrder, base, mode);
            testCase.verifySize(R,[1 numel(exponents)]);
            testCase.verifyGreaterThanOrEqual(R,0);
            testCase.verifyLessThan(R,1);
        end
        function betaBaseResidual(testCase)
            sVal = sym('sqrt(5)');
            exponents = 1:-1:-3;
            precisionOrder = 5;
            base = 2.7;
            mode = 'residual';
            R = ExpoExpand(sVal, exponents, precisionOrder, base, mode);
            testCase.verifySize(R,[1 numel(exponents)]);
            testCase.verifyGreaterThanOrEqual(R,0);
            testCase.verifyLessThan(R,1);
        end
    end
end
