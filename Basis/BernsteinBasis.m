classdef BernsteinBasis < BSplineBasis
    methods
        function B = BernsteinBasis(domain, degree)
            knots = [domain(1) * ones(1, degree + 1), domain(2) * ones(1, degree + 1)];
            B@BSplineBasis(knots, degree)
        end
    end
end