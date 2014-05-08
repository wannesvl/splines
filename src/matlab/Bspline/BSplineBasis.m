classdef BSplineBasis < Basis
    methods
        function B = BSplineBasis(knots, degree)
            % Class constructor
            B@Basis(knots, degree)
        end

        function b = f(self, x)
            % Evaluate the BSplineBasis at x
            %
            % 
            x = x(:);
            k = self.knots;
            basis = cell(self.degree + 1, 1);
            basis{1} = cell2mat(arrayfun(@(i) self.ind(i, x), ...
                                (1:length(k) - 1), ...
                                'UniformOutput', false));
            for d=1:self.degree + 1
                B = zeros(length(x), length(k) - d - 1);
                for i=1:length(k) - d - 1
                    basisd = basis{d};
                    b = 0 * x;
                    denom = k(i + d) - k(i);
                    if denom ~= 0
                        b = (x - k(i)) .* basisd(:, i) / denom;
                    end
                    denom = k(i + d + 1) - k(i + 1);
                    if denom ~= 0
                        b = b + (k(i + d + 1) - x) .* basisd(:, i + 1) / denom;
                    end
                    B(:, i) = b;
                end
                basis{d + 1} = B;
            end
            b = basis{self.degree + 1};
        end

        function [B, P] = derivative(self, o)
            if nargin == 1
                o = 1;
            end
            P = eye(length(self));
            k = self.knots;
            B = self.cl(self.knots(o:end-o), self.degree - o);
            for i=0:o-1
                k = k(2:end-1);
                delta_k = k(self.degree - i + 1:end) - k(1:end - self.degree + i);
                T = zeros(length(self) - 1 - i, length(self) - i);
                T(1:length(self) - i:end) = -1 ./ delta_k;
                T(length(self) - i:length(self) - i:end) = 1 ./ delta_k;
                P = (self.degree - i) * T * P;
            end
        end

        function T = transform(self, other)
            x = self.greville();
            if strcmp(class(other), class(self))
                T = self.f(x) \ other.f(x);
            end
            T(abs(T) < 1e-10) = 0;
        end
    end
end
