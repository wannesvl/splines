
% splines - A framework for manipulating and optimizing (B-)splines
% Copyright (C) 2015 Wannes Van Loock
%
% This file is part of splines.
%
% splines is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% splines is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU Lesser General Public License for more details.
%
% You should have received a copy of the GNU Lesser General Public License
% along with splines.  If not, see <http://www.gnu.org/licenses/>.

% classdef (Abstract) UnivariateBasis
classdef UnivariateBasis

    properties (Access=protected)
        cl  % The name of the class
    end
    properties (SetAccess={?UnivariateBasis}, GetAccess={?UnivariateBasis,?Function})
        x_  % #length(basis) of independent points
        fx_ %
    end
    % properties %(Abstract)
    %     degree
    % end

    % Each subclass should implement these methods
    % methods (Abstract)
    %     length(self)
    %     plus(self, other)
    %     mtimes(self, other)
    %     transform(self, other)
    %     f(self, x)
    % end

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
