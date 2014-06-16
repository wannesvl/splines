classdef BSplineTest < matlab.unittest.TestCase

    methods (Test)
        function testf(testCase)
            B1 = BSplineBasis([0, 0, linspace(0, 1, 11), 1, 1], 2);
            B2 = BSplineBasis([0, 0, 0, linspace(0, 1, 11), 1, 1, 1], 2);
            C = randn(length(B1), length(B2));
            T = BSpline({B1, B2}, C);
            x = linspace(0, 1, 101);
            Tf = cell2mat(T.f({x, x}));
            testCase.verifyEqual(Tf, B1.f(x) * C * B2.f(x));
        end

        
    end

end