classdef (Abstract) UnivariateBasis
    properties (Access=protected)
        cl
    end

    properties (SetAccess={?UnivariateBasis}, GetAccess={?Function})
        x_  % length of basis independent points on which the basis can be evaluated
    end

    properties (Abstract)
        degree
    end
    methods (Abstract)
        % Each subclass must define these methods
        length(self)
        plus(self, other)
        mtimes(self, other)
        transform(self, other)  % Transformation matrix from 1 basis to another
        f(self, x)  % Call the basis at x
    end
    methods
        function b = UnivariateBasis(degree)
            validateattributes(degree, {'numeric'}, {'scalar', 'integer'});
            b.degree = degree;
            b.cl = str2func(class(b));
        end
    end
end