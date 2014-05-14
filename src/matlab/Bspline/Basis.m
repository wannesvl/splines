classdef Basis
    properties (Access=protected)
        cl
    end
    properties
        knots
        degree
    end
    methods
        function B = Basis(knots, degree)
            % Constructor for Basis
            %
            % Args:
            %    knots (vector, double): the knot sequence of the basis
            %    degree (int): the degree of the basis
            %
            % Returns:
            %    Basis: an instance of the Basis class
            B.knots = knots(:);
            B.degree = degree;
            B.cl = str2func(class(B));
        end

        function s = length(self)
            s = length(self.knots) - self.degree - 1;
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
                I = (x >= self.knots(i)) .* (x <= self.knots(i+1));
            else
                I = (x > self.knots(i)) .* (x <= self.knots(i+1));
            end
        end

        function c = count_knots(self, k)
            % Count occurence of knot k
            %
            % Args:
            %    k (double): knot to count
            %
            % Returns:
            %    int: Number of occurences of the knot or nan if count == 0
            c = histc(self.knots, k);
            if c == 0
                c = nan;
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
            %    Basis: The combined basis
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
            % Returns:
            %    Basis: The sum of self and other
            if strcmp(class(other), class(self))
                degree = max(self.degree, other.degree);
                b = self.combine(other, degree);
            elseif strcmp(class(other), 'double')
                b = self;
            else
                error('Incompatible datatype')
            end
        end

        function b = mtimes(self, other)
            % Returns the product of two bases
            %
            % Returns:
            %    Basis: The product of self and other
            if strcmp(class(other), class(self))
                degree = self.degree + other.degree;
                b = self.combine(other, degree);
            elseif strcmp(class(other), 'double')
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
            g = arrayfun(@(k) sum(self.knots(k + 1:k + self.degree)) / self.degree, ...
                         (1:length(self)));
        end

        function b = insert_knots(self, knots)
            % Insert knots in the basis
            %
            % Args:
            %    knots (vector, double): The desired knot insertion locations
            %
            % Returns:
            %    Basis: The refined basis
            knots = sort([self.knots; knots(:)]);
            b = self.cl(knots, self.degree);
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
    end
end
        