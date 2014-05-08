classdef BSpline
    properties (Access=protected)
        cl
    end
    properties
        coeffs
        basis
    end
    methods
        function s = BSpline(basis, coeffs)
            % Constructor for BSpline
            s.basis = basis;
            s.coeffs = coeffs(:);
            s.cl = str2func(class(s));
        end

        function s = plus(self, other)
            if strcmp(class(self), class(other))
                basis = self.basis + other.basis;
                coeffs = basis.transform(self.basis) * self.coeffs + ...
                         basis.transform(other.basis) * other.coeffs;
            else
                basis = self.basis;
                coeffs = self.coeffs + other;
            end
            s = self.cl(basis, coeffs);
        end

        function s = uminus(self)
            s = self.cl(self.basis, -self.coeffs)
        end

        function s = minus(self, other)
            s = self + (- other)
        end

        function s = insert_knots(self, knots)
            basis = self.basis.insert_knots(knots)
            coeffs = basis.transform(self.basis) * self.coeffs
            s = self.cl(basis, coeffs)
        end

        function i = integral(self)
            k = self.basis.knots;
            d = self.basis.degree;
            i = sum((k(d + 2:end) - k(1:end - d - 1)) .* self.coeffs) / (d + 1);
        end
    end
end