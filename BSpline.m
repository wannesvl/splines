classdef BSpline < Function
    methods
        function p = BSpline(basis, coeffs)
            p@Function(basis, coeffs);
            if ~all(cellfun(@(c) isa(c, 'BSplineBasis'), p.basis))
                error('All bases are required to be polynomial basis')
            end 

        end

        function s = mtimes(self, other)
            function T = transform(A, B)
                T = A \ B;
                T(abs(T) < 1e-10) = 0;
            end
            if isa(self, class(other))
                basis = cellfun(@mtimes, self.basis, other.basis, 'UniformOutput', false);
                [i_other, i_self] = cellfun(@(b1, b2) b1.pairs(b2), self.basis, other.basis, 'UniformOutput', false);
                coeffs_product = self.coeffs(i_self{:}) * other.coeffs(i_other{:});
                % Determine transformation matrices
                x = cellfun(@(b) b.x_, basis, 'UniformOutput', false);
                b = cellfun(@(b, x) b.f(x), basis, x, 'UniformOutput', false);
                b_self = cellfun(@(b, x) b.f(x), self.basis, x, 'UniformOutput', false);
                b_other = cellfun(@(b, x) b.f(x), other.basis, x, 'UniformOutput', false);
                basis_product = cellfun(@(b1, b2, is, io) b1(:, is) .* b2(:, io), b_self, b_other, i_self, i_other, 'UniformOutput', false);
                T = cellfun(@(b, bi) transform(b, bi), b, basis_product, 'UniformOutput', false);
                s = self.cl(basis, T * coeffs_product);;
            else
                s = mtimes@Function(self, other);
            end
        end

        function s = plus(self, other)
            if isa(self, class(other))
                s = plus@Function(self, other);
            else  % Assume plus with array
                try
                    basis = self.basis;
                    coeffs = self.coeffs + other;
                    s = self.cl(basis, coeffs);
                catch err
                    s = other + self;
                end
            end
        end

        function i = integral(self)
            T = cellfun(@(b) b.integral, self.basis, 'UniformOutput', false);
            i = T * self.coeffs;
            i = i.coeffs{1};
        end
    end
end