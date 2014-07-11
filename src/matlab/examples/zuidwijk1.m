close all; clear all; clc;

d = 3;  % degree
L = 4;  % Range [0, L]
n = 21;  % Number of knots
Bl = BSplineBasis([0 * ones(1, d) linspace(0, L, n) L * ones(1, d)], d);
l = Polynomial([0, 1]);  % The parameter

% 1. Solve primal
% ===============
x = BSpline.sdpvar(Bl, [1, 2]);
obj = x(1) + 2 * x(2);
con = [x(1) >= 0, x(2) >= 0, x(2) <= 2, x(1) + l * x(2) <= 2];
options = sdpsettings('verbose',1);
sol = solvesdp(con, -obj.integral, options);

L_ = linspace(0, L, 101);
x = double(x);
obj = x(1) + 2 * x(2);
plot(L_, obj.f(L_))
hold on
plot(L_(L_<=1), 6 - 2 * L_(L_<=1), 'r')
plot(L_(L_<=2 & L_>=1), 4 ./ L_(L_<=2 & L_>=1), 'r')
plot(L_(L_>=2), 2 * ones(length(L_(L_ >= 2))), 'r')

% 2. Solve Dual
% =============
% y = BSpline.sdpvar(Bl, [1, 4]);
% obj = 2 * y(3) + 2 * y(4);
% con = [y >= 0, y(4) - y(1) == 1, -y(2) + y(3) + l * y(4) == 2];
y = BSpline.sdpvar(Bl, [1, 2]);
obj = 2 * y(1) + 2 * y(2);
con = [y >= 0, y(2) >= 1, y(1) + l * y(2) >= 2];
options = sdpsettings('verbose',1);
sol = solvesdp(con, obj.integral, options);

L_ = linspace(0, L, 101);
y = double(y);
obj = double(obj);
plot(L_, obj.f(L_), 'g')
xlabel('\theta')
ylabel('x_1(\theta) + 2 x_2(\theta)')
