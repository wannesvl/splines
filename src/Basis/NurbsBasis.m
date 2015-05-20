classdef NurbsBasis < PieceWiseBasis
    properties
        weights
    end

    methods
        function basis = NurbsBasis(knots, degree, weights)
            basis@PieceWiseBasis(knots, degree);
            basis.weights = weights(:);
        end

        function b = f(self, x)
            basis_denom = BSplineBasis(self.knots, self.degree);
            coeffs_denom = Coefficients(self.weights, size(self.weights), [1, 1]);
            denom = BSpline(basis_denom, coeffs_denom);
            b = bsxfun(@times, basis_denom.f(x), self.weights');
            b = bsxfun(@rdivide, b, denom.f(x));
        end
    end
end
