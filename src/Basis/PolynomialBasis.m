classdef PolynomialBasis < UnivariateBasis
    properties
        degree
    end
    methods
        function p = PolynomialBasis(degree)
            p@UnivariateBasis(degree);
            p.x_ = linspace(0, 1, degree + 1);
        end

        function l = length(self)
            l = self.degree + 1;
        end

        function s = plus(self, other)
            if isa(other, class(self))
                degree = max(self.degree, other.degree);
            elseif isa(other, 'numeric')  % The basis does not change
                degree = self.degree;
            elseif isa(self, 'numeric')
                degree = other.degree;
            end
            s = PolynomialBasis(degree);
        end

        function s = mtimes(self, other)
            if isa(other, class(self))
                degree = self.degree + other.degree;
            elseif isa(other, 'numeric')  % The basis does not change
                degree = self.degree;
            elseif isa(self, 'numeric')
                degree = other.degree;
            end
            s = self.cl(degree);
        end

        function T = transform(self, other)
            validateattributes(self.degree, {'numeric'}, {'>=', other.degree})
            T = zeros(length(self), length(other));
            T(1:length(other), 1:length(other)) = eye(length(other));
        end

        function b = f(self, x)
            if isa(x, 'BSpline')
                b = cell(length(self), 1);
                for i=1:length(self)
                    b{i} =  x ^ (i - 1);
                end
            else
                x = x(:);
                b = zeros(length(x), length(self));
                for i=1:length(self)
                    b(:, i) = x .^ (i - 1);
                end
            end
        end
    end
end
