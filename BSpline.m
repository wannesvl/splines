classdef BSpline < Function
    methods
        function p = BSpline(basis, coeffs)
            p@Function(basis, coeffs);
            if ~all(cellfun(@(c) isa(c, 'BSplineBasis'), p.basis))
                error('All bases are required to be BSplineBasis basis')
            end
        end

        function s = mtimes(self, other)
            function T = transform(A, B)
                T = A \ B;
                T(abs(T) < 1e-10) = 0;
            end
            if ~isa(self, mfilename)  % make sure self is the BSpline object
                temp = self;
                self = other;
                other = temp;
            end
            if isa(other, 'Polynomial')
                domain = cellfun(@(b) [b.knots(1), b.knots(end)], self.basis, 'UniformOutput', false);
                other = other.to_bspline(domain);
            end
            if isa(self, class(other))
                basis = cellfun(@mtimes, self.basis, other.basis, 'UniformOutput', false);
                [i_self, i_other] = cellfun(@(b1, b2) b1.pairs(b2), self.basis, other.basis, 'UniformOutput', false);
                coeffs_product = self.coeffs(i_self{:}) * other.coeffs(i_other{:});
                % Determine transformation matrices
                x = cellfun(@(b) b.x_, basis, 'UniformOutput', false);
                b = cellfun(@(b, x) b.f(x), basis, x, 'UniformOutput', false);
                b_self = cellfun(@(b, x) b.f(x), self.basis, x, 'UniformOutput', false);
                b_other = cellfun(@(b, x) b.f(x), other.basis, x, 'UniformOutput', false);
                basis_product = cellfun(@(b1, b2, is, io) b1(:, is) .* b2(:, io), b_self, b_other, i_self, i_other, 'UniformOutput', false);
                T = cellfun(@(b, bi) transform(b, bi), b, basis_product, 'UniformOutput', false);
                s = self.cl(basis, T * coeffs_product);
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

        function d = derivative(self, ord, coord)
            % Derivative on coord of self
            b = self.basis;
            bi = self.basis{coord};
            [dbi, P] = bi.derivative(ord);
            T = cellfun(@(p) eye(length(p)), b, 'UniformOutput', false);
            T{coord} = P;
            b{coord} = dbi;
            d = self.cl(b, T * self.coeffs);
        end

        function d = add_basis(self, basis, i)
            % Add basis to dimension i to the spline
            rep = ones(1, self.dims + 1);
            rep(i) = length(basis);
            % Repeat and permute coefficients along the new dimension
            coeffs = Coefficients(repmat(...
                                  permute(self.coeffs.coeffs, circshift((1:self.dims+1)', i)),...
                                  rep));
            % Add basis to list of bases
            b = cell(1, self.dims + 1);
            b(arrayfun(@(a) a ~= i, 1:self.dims+1)) = self.basis;
            b(i) = {basis};
            d = self.cl(b, coeffs);
        end

        function s = insert_knots(self, knots)
            % Insert knots in each basis
            %
            % Args: 
            %    knots (cell): Desired knot insertions for each
            %    basis
            %
            % Returns:            
            %    BSpline: Bspline with updated coefficients for the refined
            %    knot sequences
            b = cellfun(@(b, k) b.insert_knots(k), self.basis, knots, 'UniformOutput', false);
            T = cellfun(@(b1, b2) b1.transform(b2), b, self.basis, 'UniformOutput', false);
            s = self.cl(b, T * self.coeffs);
        end
    end
end