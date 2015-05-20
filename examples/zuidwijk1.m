
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

close all; clear all; clc;

d = 3;  % degree
L = 4;  % Range [0, L]
n = 21;  % Number of knots
Bl = BSplineBasis([0 L], d, n);
l = parameter();

% 1. Solve primal
% ===============
x = BSpline.sdpvar(Bl, [2, 1]);
obj = x(1) + 2 * x(2);
con = [x(1) >= 0, x(2) >= 0, x(2) <= 2, x(1) + l * x(2) <= 2];
options = sdpsettings('verbose',1);
sol = optimize(con, -obj.integral, options);

L_ = linspace(0, L, 101);
obj = value(obj);
plot(L_, obj.f(L_))
hold on
plot(L_(L_<=1), 6 - 2 * L_(L_<=1), 'r')
plot(L_(L_<=2 & L_>=1), 4 ./ L_(L_<=2 & L_>=1), 'r')
plot(L_(L_>=2), 2 * ones(length(L_(L_ >= 2))), 'r')

% 2. Solve Dual
% =============
y = BSpline.sdpvar(Bl, [1, 2]);
obj = 2 * y(1) + 2 * y(2);
con = [y >= 0, y(2) >= 1, y(1) + l * y(2) >= 2];
options = sdpsettings('verbose',1);
sol = optimize(con, obj.integral, options);

L_ = linspace(0, L, 101);
y = value(y);
obj = value(obj);
plot(L_, obj.f(L_), 'g')
xlabel('\theta')
ylabel('x_1(\theta) + 2 x_2(\theta)')
