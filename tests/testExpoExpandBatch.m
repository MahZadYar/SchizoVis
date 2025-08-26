classdef testExpoExpandBatch < matlab.unittest.TestCase
    methods(Test)
        function testSmallDigitsExtraction(testCase)
            n = 5; base = 10; precisionOrder = 5; mode = "digits"; [~, sVal] = SchizoGen(n, base);
            p_min = -precisionOrder; exponents = 2:-1:p_min; % simple small window
            D = expoExpandBatch(sVal, exponents, precisionOrder, base, mode);
            testCase.verifySize(D, [1 numel(exponents)]);
            testCase.verifyTrue(all(D>=0 & D < base));
        end
            function basicDigits(testCase)
                sVal = sym('sqrt(2)');
                exponents = [1 0 -1 -2];
                precisionOrder = 3;
                base = 10;
                mode = 'digits';
                D = expoExpandBatch(sVal, exponents, precisionOrder, base, mode);
                testCase.verifySize(D,[1 numel(exponents)]);
                testCase.verifyGreaterThanOrEqual(D,0);
            end
end
