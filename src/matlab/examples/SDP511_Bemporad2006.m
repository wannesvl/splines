close all; clear all; clc;

% Problem data
% ============
A0 = [1 2 -3; 2 4 -1; -3 -1 3];
At1 = [1 -1 2; -1 1 3; 2 3 2];
At2 = [-1 1 0; 1 1 2; 0 2 -2];
Ax1 = [3 -2 4; -2 1 -2; 4 -2 -2];
Ax2 = [-3 1 1; 1 -2 -1; 1 -1 1];
Ax3 = [5 4 2; 4 1 1; 2 1 -1];

% 1. Determine Feasible region
% ============================
degree = 3;
knots = [0 * ones(1, degree), linspace(0, 2 * pi, 20), 2 * pi * ones(1, degree)];
step = 1 / 10;
knots = 2 * pi * ((-degree * step):step:(1 + degree * step));
b_phi = BSplineBasis(knots, degree);
phi = linspace(0, 2 * pi, 501)';

% Optimization problem
t = BSpline.sdpvar(b_phi, [1, 2]);
x = BSpline.sdpvar(b_phi, [1, 3]);
K = A0 + At1 * t(1) + At2 * t(2) + Ax1 * x(1) + Ax2 * x(2) + Ax3 * x(3);
options = sdpsettings('verbose', 1, 'solver', 'sdpt3');
con = [K >= 0, t(2) >= -2, t(1) <= 2];
for i=1:degree  % Continuity constraints
    con = [con, t.coeffs.coeffs{i} == t.coeffs.coeffs{end-degree+i}];
end
sol = solvesdp(con, -sum(t(1).f(phi) .* sin(phi) + t(2).f(phi) .* cos(phi)), options);

% Plot feasible region
t = double(t);
figure
phi = linspace(0, 2 * pi, 501);
plot(t(1).f(phi), t(2).f(phi))

% 2. Solve optimization problem
% =============================
% fill boundary
% dirty hacking... Maybe we can automate this?
b_beta = BSplineBasis([0, 0, 1, 1], 1);
origin = [0, 0];
c = cell(length(b_phi), 1);
[c{:}] = deal(origin);
ct = t.coeffs - origin;
t = BSpline({b_phi, b_beta}, Coefficients([c, ct.coeffs]))

% Define optimization problem
degree = 3;
na = 5;
nb = 10;
Ba = b_phi;
Bb = BSplineBasis([0 * ones(1, degree), linspace(0, 1, nb), ones(1, degree)], degree);

x = BSpline.sdpvar({Ba, Bb}, [3, 1]);
K = A0 + At1 * t(1) + At2 * t(2) + Ax1 * x(1) + Ax2 * x(2) + Ax3 * x(3);
objective = x(1) - 2 * x(2) + x(3);
sol = solvesdp([K >= 0], objective.integral, options);

% Plotting
a = linspace(0, 2 * pi, 101);
b = linspace(0, 1, 101);
T1 = t(1).f({a, b});
T2 = t(2).f({a, b});
x = double(x);
objective = x(1) - 2 * x(2) + x(3);
O_approx = objective.f({a, b});
figure
surf(T1, T2, O_approx, 'EdgeColor', 'none')
camlight left; light; lighting phong; shading interp; alpha(0.5);
return
% Solve gridded problem
% =====================
O = zeros(21, 21);
x = sdpvar(3, 1);
t = sdpvar(2, 1);
K = A0 + At1 * t(1) + At2 * t(2) + Ax1 * x(1) + Ax2 * x(2) + Ax3 * x(3);
objective = x(1) - 2 * x(2) + x(3);
prob = optimizer([K >= 0], objective, options, t, x);
for i=1:5:length(a)
    for j=1:5:length(b)
        [xk, errorcode] = prob{[T1(i, j); T2(i, j)]};
        O(ceil(i / 5), ceil(j / 5)) = xk(1) - 2 * xk(2) + xk(3);
    end
end

% Plot results
hold on;
surf(T1(1:5:length(a), 1:5:length(b)), T2(1:5:length(a), 1:5:length(b)), O)
camlight left; light; lighting phong; shading interp; alpha(0.5);

figure
surf(T1(1:5:length(a), 1:5:length(b)), T2(1:5:length(a), 1:5:length(b)), O - O_approx(1:5:length(a), 1:5:length(b)))
