classdef BSplineBasis < PieceWiseBasis
    methods
        function basis = BSplineBasis(knots, degree, n)
            % Constructor for BSplineBasis, a subclass of PieceWiseBasis
            %
            % Args:
            %    knots (vector, double): the knot sequence of the basis
            %    degree (int): the degree of the basis
            %    n (int): optional argument for easy equidistant knot choice
            %
            % Returns:
            %    BSplineBasis: an instance of the Basis class
            %
            % Example:
            %    > B = BsplineBasis([0, 0, 0, 0.5, 1, 1, 1], 2)
            if nargin == 3
                knots = [knots(1) * ones(1, degree), ...
                         linspace(knots(1), knots(2), n), ...
                         knots(2) * ones(1, degree)];
            end
            basis@PieceWiseBasis(knots, degree)
        end

        function b = f(self, x)
            % Evaluate the BSplineBasis at x
            %
            % This function implements the Cox-de Boor recursive formula for
            % the evaluation of B-spline basis functions
            %
            % Args:
            %    x (vector, double): The evaluation sites
            %
            % Returns:
            %    matrix, double: Each column of b contains the evaluated basis
            %        function at the evaluation sites
            x = x(:);
            k = self.knots;
            % basis = cell(self.degree + 1, 1);
            % basis{1} = cell2mat(arrayfun(@(i) self.ind(i, x), ...
            %                     (1:length(k) - 1), ...
            %                     'UniformOutput', false));
            B = zeros(length(x), length(k) - 1);
            for i=1:length(k) - 1
                B(:, i) = self.ind(i, x);
            end
            basisd = B;
            sx = size(x);
            for d=1:self.degree
                B = zeros(length(x), length(k) - d - 1);
                for i=1:length(k) - d - 1
                    b = zeros(sx);
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
                basisd = B;
            end
            b = sparse(B);
        end

        function [B, P] = derivative(self, o)
            if nargin == 1
                o = 1;
            end
            P = eye(length(self));
            k = self.knots;
            B = self.cl(self.knots(o + 1:end - o), self.degree - o);
            for i=0:o-1
                k = k(2:end-1);
                delta_k = k(self.degree - i + 1:end) - k(1:end - self.degree + i);
                T = zeros(length(self) - 1 - i, length(self) - i);
                T(1:length(self) - i:end) = -1 ./ delta_k;
                T(length(self) - i:length(self) - i:end) = 1 ./ delta_k;
                P = (self.degree - i) * T * P;
            end
            P = sparse(P);
        end

        function i = integral(self)
        % Integral of the basis functions
            k = self.knots;
            d = self.degree;
            i = (k(d + 2:end) - k(1:end - d - 1))'  / (d + 1);
        end
    end
end
