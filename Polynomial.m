classdef Polynomial < Function
    methods
        function p = Polynomial(varargin)
            if nargin == 1  % Assume only the coefficients are given
                coeffs = varargin{1};
                if ~isa(coeffs, 'cell')  % Scalar valued polynomial
                    coeffs = num2cell(coeffs);
                end
                if isvector(coeffs)
                    coeffs = coeffs(:);
                    basis = {PolynomialBasis(length(coeffs) - 1)};
                else  % And assume coeffs is correctly formatted
                    basis = arrayfun(@PolynomialBasis, size(coeffs) - 1, 'UniformOutput', false);
                end
            else
                basis = varargin{1};
                coeffs = varargin{2};
            end
            p@Function(basis, coeffs);
            if ~all(cellfun(@(c) isa(c, 'PolynomialBasis'), p.basis))
                error('All bases are required to be polynomial basis')
            end
        end

        function s = plus(self, other)
            if isa(self, class(other))
                s = plus@Function(self, other);
            else  % Assume plus with array
                try
                    basis = self.basis; 
                    coeffs = self.coeffs;
                    coeffs.coeffs{1} = coeffs.coeffs{1} + other;
                    s = self.cl(basis, coeffs);
                catch err
                    s = other + self;
                end
            end
        end

        function s = to_bspline(self, domain)
            % Convert the polynomial to a BSpline representation on the
            % rectangular domain
            %
            % TODO: Add polygonal domains
            convert = @(b, d) BSplineBasis([d(1) * ones(1, b.degree + 1), ...
                                            d(end) * ones(1, b.degree + 1)], ...
                                            b.degree);
            basis = cellfun(convert, self.basis, domain, 'UniformOutput', false);
            T = cellfun(@(b1, b2) b1.transform(b2), basis, self.basis, 'UniformOutput', false);
            s = BSpline(basis, T * self.coeffs);
        end
    end
end