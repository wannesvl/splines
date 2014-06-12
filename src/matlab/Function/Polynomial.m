classdef Polynomial < Function
    methods
        function p = Polynomial(basis, coeffs)
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
    end
end