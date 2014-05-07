classdef Basis
    properties (Access=protected)
        cl
    properties
        knots
        degree
    end
    methods
        function B = Basis(knots, degree)
            % Constructor
            B.knots = knots;
            B.degree = degree;
            B.cl = str2func(class(B));
        end

        function s = size(self)
            s = size(self.knots) - self.degree - 1;
        end

        function i = ind(self, i, x)
            % indicator function between knots(i) and knots(i+1)
            if i < self.degree + 1 && self.knots(1) == self.knots(i)
                i = (x >= self.knots(i)) .* (x <= self.knots(i+1));
            else
                i = (x > self.knots(i)) .* (x <= self.knots(i+1));
            end
        end

        function c = count_knots(self)
            % count occurence of each knot
            k = unique(self.knots);
            c = containers.Map(k, histc(self.knots, k));
        end

        function b = combine(self, other, degree)
            % Combine 2 bases to form a new basis of degree 
            c_self = self.count_knots;
            c_other = other.count_knots;
            breaks = union(self.knots, other.knots);
            knots = [];
            for i = 1:length(breaks)
                try
                    mi = c_self(breaks(i));
                catch err
                    mi = 0;
                end
                try
                    ni = c_other(breaks(i));
                catch err
                    ni = 0;
                end
                mult = max(mi + degree - self.degree, ni + degree - other.degree);
                knots = [knots; ones(mult, 1) * breaks(i)];
            end
            cl = str2func(class(self));
            b = cl(knots, degree);
        end

        function b = plus(self, other)
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
            g = 0;
        end

        function b = insert_knots(self, knots)
            knots = sort([self.knots, knots]);
            cl = str2func(class(self));
            b = cl(knots, self.degree);
        end
    end
end
        