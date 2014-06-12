classdef BernsteinBasis < BSplineBasis
    methods
        function B = BernsteinBasis(domain, degree)
            knots = [domain(1) * ones(1, degree), ...
                     linspace(domain(1), domain(2)), ...
                     domain(2) * ones(1, degree)];
            B@BSplineBasis(knots, degree)
        end
    end
end