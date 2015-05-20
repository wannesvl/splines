
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

function p = parameter(n)
    % Creates n linear polynomial functions in n variables
    if nargin == 0 || n == 1
        p = Polynomial({0, 1}');
        return
    end
    data = zeros(n, 2^n);
    data(:, 1 + 2.^(0:n-1)) = eye(n);
    coeffs = Coefficients(data, [1, 2^n], [n, 1]);
    I = reshape(1:2^n, 2 * ones(1, n));
    coeffs = coeffs(I);
    p = Polynomial(coeffs);
