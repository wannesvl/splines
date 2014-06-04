classdef TensorSpline
    properties (Access=protected)
        cl
    end
    properties
        coeffs
        basis
    end
    methods
        function s = TensorSpline(basis, coeffs)
            % Constructor for TensorSpline
            %
            % Only supports scalar coefficients for the moment
            s.basis = basis;
            if isa(coeffs, 'BSplineCoeffs')
                s.coeffs = coeffs;
            else
                s.coeffs = BSplineCoeffs(coeffs);
            end
            % Validate input
            
            % if length(s.basis) ~= length(s.coeffs)
            %     error('B-spline coefficient of different size than basis')
            % end
            s.cl = str2func(class(s));
        end

        function s = f(self, x)
            s = cellfun(@(b, x) b.f(x), self.basis, x, 'UniformOutput', false) * self.coeffs;
        end

        function d = dims(self)
            d = length(self.basis);
        end

        function s = plus(self, other)
            if isa(other, 'TensorSpline')
                basis = cellfun(@plus, self.basis, other.basis, 'UniformOutput', false);
                Tself = cellfun(@(b1, b2) b1.transform(b2), basis, self.basis, 'UniformOutput', false);
                Tother = cellfun(@(b1, b2) b1.transform(b2), basis, other.basis, 'UniformOutput', false);
                coeffs = Tself * self.coeffs + Tother * other.coeffs;
            else
                basis = self.basis;
                coeffs = self.coeffs + other;
            end
            s = self.cl(basis, coeffs);
        end
    end
end
