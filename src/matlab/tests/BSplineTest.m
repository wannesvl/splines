classdef BSplineTest < matlab.unittest.TestCase

    methods (Test)
        function testf(testCase)
            B1 = BSplineBasis([0, 0, linspace(0, 1, 11), 1, 1], 2);
            B2 = BSplineBasis([0, 0, 0, linspace(0, 1, 11), 1, 1, 1], 3);
            C = randn(length(B1), length(B2));
            T = BSpline({B1, B2}, C);
            x = linspace(0, 1, 101);
            Tf = T.f({x, x});
            testCase.verifyEqual(Tf, B1.f(x) * C * B2.f(x)', 'AbsTol', 1e-14);
        end

        function testmtimes(testCase)
            B1 = BSplineBasis([0, 0, linspace(0, 1, 11), 1, 1], 2);
            B2 = BSplineBasis([0, 0, 0, linspace(0, 1, 11), 1, 1, 1], 3);
            C1 = randn(length(B1), length(B2));
            C2 = randn(length(B1), length(B2));
            T1 = BSpline({B1, B2}, C1);
            T2 = BSpline({B1, B2}, C2);
            x = linspace(0, 1, 101);
            T1f = T1.f({x, x});
            T2f = T2.f({x, x});
            G = T1 * T2;
            testCase.verifyEqual(T1f .* T2f, G.f({x, x}), 'AbsTol', 1e-14);
        end

        function testplus(testCase)
            B1 = BSplineBasis([0, 0, linspace(0, 1, 11), 1, 1], 2);
            B2 = BSplineBasis([0, 0, 0, linspace(0, 1, 11), 1, 1, 1], 3);
            C1 = randn(length(B1), length(B2));
            C2 = randn(length(B1), length(B2));
            T1 = BSpline({B1, B2}, C1);
            T2 = BSpline({B1, B2}, C2);
            x = linspace(0, 1, 101);
            T1f = T1.f({x, x});
            T2f = T2.f({x, x});
            G = T1 + T2;
            testCase.verifyEqual(T1f + T2f, G.f({x, x}), 'AbsTol', 1e-14);
        end
    end
end