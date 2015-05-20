
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

classdef NurbsBasis < PieceWiseBasis
    properties
        weights
    end

    methods
        function basis = NurbsBasis(knots, degree, weights)
            basis@PieceWiseBasis(knots, degree);
            basis.weights = weights(:);
        end

        function b = f(self, x)
            basis_denom = BSplineBasis(self.knots, self.degree);
            coeffs_denom = Coefficients(self.weights, size(self.weights), [1, 1]);
            denom = BSpline(basis_denom, coeffs_denom);
            b = bsxfun(@times, basis_denom.f(x), self.weights');
            b = bsxfun(@rdivide, b, denom.f(x));
        end
    end
end
