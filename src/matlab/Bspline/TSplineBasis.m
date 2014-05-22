classdef TSplineBasis < Basis
    methods
        function T = TSplineBasis(knots, degree)
            % Constructor for TSplineBasis, a subclass of Basis
            %
            % Args:
            %    knots (vector, double): the knot sequence of the basis
            %    degree (int): the degree of the basis
            %
            % Returns:
            %    TSplineBasis: an instance of the TSplineBasis class
            %
            % Example:
            %    > B = BsplineBasis([0, 0, 0, 0.5, 1, 1, 1], 2)
            T@Basis(knots, degree)
        end

        function b = f(self, x)
            % Evaluate the TSplineBasis at x
            %
            % This function implements the recursive formula for
            % the evaluation of T-spline basis functions
            %
            % Args:
            %    x (vector, double): The evaluation sites
            %
            % Returns:            
            %    matrix, double: Each column of b contains the evaluated basis
            %        function at the evaluation sites
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
                    denom = sin(0.5 * (k(i + d) - k(i)));
                    if denom ~= 0
                        b = sin(0.5 * (x - k(i))) .* basisd(:, i) / denom;
                    end
                    denom = sin(0.5 * (k(i + d + 1) - k(i + 1)));
                    if denom ~= 0
                        b = b + sin(0.5 * (k(i + d + 1) - x)) .* basisd(:, i + 1) / denom;
                    end
                    b(isnan(b)) = 0;
                    B(:, i) = b;
                end
                basis{d + 1} = B;
            end
            b = basis{self.degree + 1};
        end
    end
end
