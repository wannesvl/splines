classdef Polynomial < Function
    methods
        function p = Polynomial(varargin)
            if nargin == 1  % Assume only the coefficients are given
                coeffs = varargin{1};
                if ~isa(coeffs, 'Coefficients') && ~isa(coeffs, 'cell')  % Scalar valued polynomial
                    coeffs = num2cell(coeffs);
                end
                if isa(coeffs, 'Coefficients')
                    basis = arrayfun(@PolynomialBasis, size(coeffs) - 1, 'UniformOutput', false);
                elseif isvector(coeffs)
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
            if isa(self, 'Function') && isa(other, 'Function')
                if ~isa(other, mfilename)
                    s = other + self;  % Let the other class handle summation
                    return
                end
                s = plus@Function(self, other);
            else  % Assume plus with array
                try
                    basis = self.basis;
                    coeffs = self.coeffs;
                    coeffs(1) = coeffs(1).data + other;
                    s = self.cl(basis, coeffs);
                catch err
                    s = other + self;
                end
            end
        end

        function s = mtimes(self, other)
            % This is a general implementation but could be performed more efficiently in subclasses
            % It is recommended to overload mtimes in the subclass
            function T = transform(A, B)
                lA = length(A);
                lB = length(B);
                T = zeros(lA + lB - 1, lA * lB);
                for i = 1:lA
                    T(i:lB-1+i, lB*(i-1)+1:lB*i) = eye(lB);
                end
            end
            if isa(self, class(other))
                % Basis of product
                basis = cellfun(@mtimes, self.basis, other.basis, 'UniformOutput', false);
                % Take kronecker product of coefficients
                [i_other, i_self] = arrayfun(@(i) find(ones(size(other.coeffs, i), size(self.coeffs, i))), 1:self.dims, 'UniformOutput', false);  % Give all indices of products
                coeffs_product = self.coeffs(i_self{:}) * other.coeffs(i_other{:});
                % self.coeffs.data, i_self{:}
                % Determine transformation matrices
                T = cellfun(@(b1, b2) transform(b1, b2), self.basis, other.basis, 'UniformOutput', false);
                s = self.cl(basis, T * coeffs_product);
            elseif isa(other, 'BSpline')
                s = (other' * self')';
            else       % Assume multiplication with array
                s = mtimes@Function(self, other);
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
