close all
clear all
clc

double_ = @(c) cellfun(@double, c, 'uni', false);  % Shortcut for taking double of cell array
d = 3;  % degree
L = 4;  % Range [0, L]
n = 21;  % Number of knots
Bl = BSplineBasis([0 * ones(1, d) linspace(0, L, n) L * ones(1, d)], d);
l = Polynomial([0, 1]);  % The parameter

% 1. Solve primal
% ===============
c = sdpvar(2 * ones(1, length(Bl)), ones(1, length(Bl)));
x = BSpline(Bl, c);
obj = x(1) + 2 * x(2);
con = [x(1) >= 0, x(2) >= 0, x(2) <= 2, x(1) + l * x(2) <= 2];
options = sdpsettings('verbose',1);
sol = solvesdp(con, -obj.integral, options);

L_ = linspace(0, L, 101);
x = BSpline(Bl, double_(c));
obj = x(1) + 2 * x(2);
plot(L_, obj.f(L_))
hold on
plot(L_(L_<=1), 6 - 2 * L_(L_<=1), 'r')
plot(L_(L_<=2 & L_>=1), 4 ./ L_(L_<=2 & L_>=1), 'r')
plot(L_(L_>=2), 2 * ones(length(L_(L_ >= 2))), 'r')

% 2. Solve Dual
% =============
c = sdpvar(4 * ones(1, length(Bl)), ones(1, length(Bl)));
y = BSpline(Bl, c);
obj = 2 * y(3) + 2 * y(4);
con = [y >= 0, y(4) - y(1) == 1, -y(2) + y(3) + l * y(4) == 2];
options = sdpsettings('verbose',1);
sol = solvesdp(con, obj.integral, options);

L_ = linspace(0, L, 101);
y = BSpline(Bl, double_(c));
obj = 2 * y(3) + 2 * y(4);
plot(L_, obj.f(L_), 'g')
xlabel('\theta')
ylabel('x_1(\theta) + 2 x_2(\theta)')
