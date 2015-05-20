
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

classdef PolynomialBasis < UnivariateBasis
    properties
        degree
    end
    methods
        function p = PolynomialBasis(degree)
            p@UnivariateBasis(degree);
            p.x_ = linspace(0, 1, degree + 1);
        end

        function l = length(self)
            l = self.degree + 1;
        end

        function s = plus(self, other)
            if isa(other, class(self))
                degree = max(self.degree, other.degree);
            elseif isa(other, 'numeric')  % The basis does not change
                degree = self.degree;
            elseif isa(self, 'numeric')
                degree = other.degree;
            end
            s = PolynomialBasis(degree);
        end

        function s = mtimes(self, other)
            if isa(other, class(self))
                degree = self.degree + other.degree;
            elseif isa(other, 'numeric')  % The basis does not change
                degree = self.degree;
            elseif isa(self, 'numeric')
                degree = other.degree;
            end
            s = self.cl(degree);
        end

        function T = transform(self, other)
            validateattributes(self.degree, {'numeric'}, {'>=', other.degree})
            T = zeros(length(self), length(other));
            T(1:length(other), 1:length(other)) = eye(length(other));
        end

        function b = f(self, x)
            if isa(x, 'BSpline')
                b = cell(length(self), 1);
                for i=1:length(self)
                    b{i} =  x ^ (i - 1);
                end
            else
                x = x(:);
                b = zeros(length(x), length(self));
                for i=1:length(self)
                    b(:, i) = x .^ (i - 1);
                end
            end
        end
    end
end
