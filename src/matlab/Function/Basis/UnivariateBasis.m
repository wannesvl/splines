classdef (Abstract) UnivariateBasis

    properties (Access=protected)
        cl  % The name of the class
    end
    properties (SetAccess={?UnivariateBasis}, GetAccess={?UnivariateBasis,?Function})
        x_  % #length(basis) of independent points 
    end
    properties (Abstract)
        degree
    end

    % Each subclass should implement these methods
    methods (Abstract)
        length(self)
        plus(self, other)
        mtimes(self, other)
        transform(self, other)
        f(self, x)
    end

    methods
        function basis = UnivariateBasis(degree)
            % An abstract base class for a general function basis.
            % A basis is defined by its integer valued degree
            % 
            % Each subclass must overload the following methods:
            %   * length: return the number of basis functions
            %   * plus: return a new basis formed by the sum of two bases 
            %   * mtimes: return a new basis formed by the product of two bases 
            %   * transform: return a transformation matrix from one basis to another
            %   * f: return the evaluation of the basis at points x
            validateattributes(degree, {'numeric'}, {'scalar', 'integer'});
            basis.degree = degree;
            basis.cl = str2func(class(basis));
        end
    end
end