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

        function s = f(self, x)
            % Evaluate B-spline at x
            s = self.basis.f(x) * self.coeffs;
        end

        function s = plus(self, other)
            if strcmp(class(self), class(other))
                basis = self.basis + other.basis;
                coeffs = basis.transform(self.basis) * self.coeffs + ...
                         basis.transform(other.basis) * other.coeffs;
                s = self.cl(basis, coeffs);
            else
                try
                    basis = self.basis;
                    coeffs = other + self.coeffs;
                    s = self.cl(basis, coeffs);
                catch err
                    basis = other.basis;
                    coeffs = self + other.coeffs;
                    s = other.cl(basis, coeffs);
                end
            end
        end

        function s = uminus(self)
            s = self.cl(self.basis, -self.coeffs);
        end

        function s = minus(self, other)
            s = self + (- other);
        end

        function s = mtimes(self, other)
            if strcmp(class(self), class(other))
                basis = self.basis * other.basis;
                grev = basis.greville;
                b_self = self.basis.f(grev);
                b_other = other.basis.f(grev);
                [i, j] = self.basis.pairs(other.basis);
                basis_product = b_self(:, i) .* ...
                                b_other(:, j);
                T = basis.f(grev) \ basis_product;
                T(abs(T) < 1e-10) = 0;
                coeffs_product = self.coeffs(i) .* other.coeffs(j);
                s = self.cl(basis, T * coeffs_product);
            else
                try
                    basis = self.basis;
                    coeffs = other * self.coeffs;
                    s = self.cl(basis, coeffs);
                catch err
                    basis = other.basis;
                    coeffs = self * other.coeffs;
                    s = other.cl(basis, coeffs);
                end
            end
            
        end

        function d = derivative(self, o)
            if nargin == 1
                o = 1;
            end
            B, P = self.basis.derivative(o);
            d = self.cl(B, P * self.coeffs);
        end

        function s = insert_knots(self, knots)
            basis = self.basis.insert_knots(knots);
            coeffs = basis.transform(self.basis) * self.coeffs;
            s = self.cl(basis, coeffs);
        end

        function i = integral(self)
            k = self.basis.knots;
            d = self.basis.degree;
            i = sum((k(d + 2:end) - k(1:end - d - 1)) .* self.coeffs) / (d + 1);
        end
    end
end