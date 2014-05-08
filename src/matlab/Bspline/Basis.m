classdef Basis
    properties (Access=protected)
        cl
    properties
        knots
        degree
    end
    methods
        function B = Basis(knots, degree)
            % Constructor for Basis
            %
            % Args:
            %    knots (): the knot sequence of the basis
            %    degree (int): the degree of the basis
            %
            % Returns:
            %    Basis: an instance of the Basis class
            B.knots = knots(:);
            B.degree = degree;
            B.cl = str2func(class(B));
        end

        function s = size(self)
            s = size(self.knots) - self.degree - 1;
        end

        function I = ind(self, i, x)
            % Indicator function between knots(i) and knots(i+1)
            %
            % Args:
            %    i (int): the index of the indicator function
            %    x (): the evaluation sites of the indicator function
            %
            % Returns:
            %    : 1 if knots(i) < x <= knots(i + 1) else 0
            if i < self.degree + 1 && self.knots(1) == self.knots(i)
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
            %    int: Number of occurences of the knot
            c = histc(self.knots, k)
            % k = unique(self.knots);
            % c = containers.Map(k, histc(self.knots, k));
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

            % c_self = self.count_knots;
            % c_other = other.count_knots;
            breaks = union(self.knots, other.knots);
            knots = [];
            for i = 1:length(breaks)
                % try
                %     mi = c_self(breaks(i));
                % catch err
                %     mi = 0;
                % end
                % try
                %     ni = c_other(breaks(i));
                % catch err
                %     ni = 0;
                % end
                mi = self.count_knots(breaks(i))
                ni = other.count_knots(breaks(i))
                mult = max(mi + degree - self.degree, ni + degree - other.degree);
                knots = [knots; ones(mult, 1) * breaks(i)];
            end
            cl = str2func(class(self));
            b = cl(knots, degree);
        end

        function b = plus(self, other)
            % Returns the sum of two bases
            %
            % Returns:
            %    Basis
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
            %    Basis
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
            %    
            g = 0;
        end

        function b = insert_knots(self, knots)
            % Insert knots in the basis
            %
            % Args:
            %    knots (): The desired knot insertion locations
            %
            % Returns:
            %    Basis: The refined basis
            knots = sort([self.knots; knots(:)]);
            cl = str2func(class(self));
            b = cl(knots, self.degree);
        end
    end
end
        