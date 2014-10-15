classdef BernsteinBasis < BSplineBasis
    methods
        function basis = BernsteinBasis(domain, degree)
            knots = [domain(1) * ones(1, degree + 1), domain(2) * ones(1, degree + 1)];
            % basis@PieceWiseBasis(knots, degree)
            basis@BSplineBasis(knots, degree)
        end
    end
end