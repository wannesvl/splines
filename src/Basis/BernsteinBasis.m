
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

classdef BernsteinBasis < BSplineBasis
    methods
        function basis = BernsteinBasis(domain, degree)
            knots = [domain(1) * ones(1, degree + 1), domain(2) * ones(1, degree + 1)];
            % basis@PieceWiseBasis(knots, degree)
            basis@BSplineBasis(knots, degree)
        end
    end
end