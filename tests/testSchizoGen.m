classdef testSchizoGen < matlab.unittest.TestCase
    methods(Test)
        function testBase10Simple(testCase)
            [fVals, ~] = SchizoGen(5, 10);
            testCase.verifyEqual(double(fVals), 12345);
        end
        function testVectorInputOrder(testCase)
            nVec = [3 1 4]; [fVals, ~] = SchizoGen(nVec, 10);
            testCase.verifyEqual(double(fVals(1)), 123);
            testCase.verifyEqual(double(fVals(2)), 1);
            testCase.verifyEqual(double(fVals(3)), 1234);
        end
    end
end
