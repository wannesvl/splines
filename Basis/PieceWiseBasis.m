classdef PieceWiseBasis < UnivariateBasis
    properties
        knots
        degree
    end
    methods
        function basis = PieceWiseBasis(knots, degree)
            % An abstract base class for univariate piecewise basis functions.
            % A piecewise basis is defined by a nondecreasing sequence of
            % knots and a degree.
            %
            % Args:
            %    knots (vector, double): the (nondecreasing) knot sequence of the basis
            %    degree (int): the degree of the basis
            %
            % Returns:
            %    PieceWiseBasis: an instance of the Basis class
            basis@UnivariateBasis(degree);
            validateattributes(knots, {'numeric'}, {'nondecreasing'})
            basis.knots = knots(:);
            basis.x_ = basis.greville;
        end

        function l = length(self)
            l = length(self.knots) - self.degree - 1;
        end

        function I = ind(self, i, x)
            % Indicator function between knots(i) and knots(i+1)
            %
            % Args:
            %    i (int): the index of the indicator function
            %    x (vector, double): the evaluation sites of the indicator function
            %
            % Returns:
            %    vector, boolean: 1 if knots(i) < x <= knots(i + 1) else 0
            x = x(:);
            if i < self.degree + 2 && self.knots(1) == self.knots(i)
                % To avoid difficulties with the first knot
                I = (x >= self.knots(i)) .* (x <= self.knots(i+1));
            else
                I = (x > self.knots(i)) .* (x <= self.knots(i+1));
            end
        end

        function count = count_knots(self, k)
            % Count occurence of knot k
            %
            % Args:
            %    k (double): knot to count
            %
            % Returns:
            %    int: Number of occurences of the knot or nan if count == 0
            count = histc(self.knots, k);
            if count == 0
                count = nan;
            end
        end

        function b = combine(self, other, degree)
            % Combine 2 bases to form a new basis of degree
            %
            % Args:
            %    other (Basis): other basis to combine
            %    degree (int): Desired degree of the combination
            %
            % Returns:
            %    PieceWiseBasis: The combined basis
            breaks = union(self.knots, other.knots);
            knots = [];
            for i = 1:length(breaks)
                mi = self.count_knots(breaks(i));
                ni = other.count_knots(breaks(i));
                mult = max(mi + degree - self.degree, ni + degree - other.degree);
                knots = [knots; ones(mult, 1) * breaks(i)];
            end
            b = self.cl(knots, degree);
        end

        function b = plus(self, other)
            % Returns the sum of two bases
            %
            % Args:
            %    other (Basis)
            %
            % Returns:
            %    PieceWiseBasis: The sum of self and other
            if isa(other, class(self))
                degree = max(self.degree, other.degree);
                b = self.combine(other, degree);
            elseif isa(other, 'double')
                b = self;
            elseif isa(self, 'double')                
                b = other;
            else
                error('Incompatible datatype')
            end
        end

        function b = mtimes(self, other)
            % Returns the product of two bases
            %
            % Args:
            %    other (Basis)
            %
            % Returns:
            %    PieceWiseBasis: The product of self and other
            if isa(other, class(self)) || isa(self, class(other))
                degree = self.degree + other.degree;
                b = self.combine(other, degree);
            elseif isa(other, 'double')
                b = self;
            else
                error('Incompatible datatype')
            end
        end

        function g = greville(self)
            % Return the greville abscissae of the basis
            %
            % Returns:
            %    vector, double: The greville abscissae
            d = self.degree;
            if d == 0
                d = 1;
            end
            g = arrayfun(@(k) sum(self.knots(k + 1:k + d)) / d, ...
                         (1:length(self)));
        end

        function b = insert_knots(self, knots)
            % Insert knots in the basis
            %
            % Args:
            %    knots (vector, double): The desired knot insertion locations
            %
            % Returns:
            %    PieceWiseBasis: The refined basis
            knots = sort([self.knots; knots(:)]);
            b = self.cl(knots, self.degree);
        end

        function basis = increase_degree(self, d)
            % Increase the degree of the basis by d
            %
            % Args:
            %    d (int): The desired degree increase
            %
            % Returns:
            %    PieceWiseBasis: The new basis
            degree = self.degree + d;
            knots = sort([self.knots; repmat(unique(self.knots), d, 1)]);
            basis = self.cl(knots, degree);
        end

        function s = support(self)
            % Return a matrix of support intervals for each basis function
            %
            % Returns:
            %    (n x 2) matrix, double: the support intervals
            s = [self.knots(1:end - self.degree - 1), ...
                 self.knots(self.degree + 2:end)];
        end

        function [i, j] = pairs(self, other)
            % Returns indices of nonzero products of basis functions
            %
            % Returns:
            %    vector, int: the valid indices of self
            %    vector, int: the corresponding valid indices of other
            is_valid = @(a, b) max(a(1), b(1)) < min(a(2), b(2));
            s_self = self.support;
            s_other = other.support;
            p = zeros(length(self), length(other));
            for i = 1:length(self)
                for j = 1:length(other)
                    p(i, j) = is_valid(s_self(i, :), s_other(j, :));
                end
            end
            [i, j, dummy] = find(p);
        end

        function T = transform(self, other)
            % Returns a transformation matrix T from self to other such that
            % self * T = other. The transformation matrix is only defined if
            % other is included in self.
            %
            % Args:
            %    other (PieceWiseBasis): Transforming basis.
            %
            % Returns:
            %    array: the transformation matrix
            x = self.x_;
            if any(diff(x) == 0)  % Fix bad x
                x = linspace(x(1), x(end), 501);
            end
            T = self.f(x) \ other.f(x);
            if any(isnan(T))
                error('Transformation matrix cannot be determined')
            end
            T(abs(T) < 1e-10) = 0;
            T = sparse(T);
        end
    end
end
        