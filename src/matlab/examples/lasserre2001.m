close all; clear all; clc;
options = sdpsettings('verbose', 1, 'solver', 'sdpt3');
% Polynomial optimization Example Lasserre2001

% bases
degree = 4;
n = 8;
B1 = BSplineBasis([0 * ones(1, degree), linspace(0, 3, n), 3 * ones(1, degree)], degree);
B2 = BSplineBasis([0 * ones(1, degree), linspace(0, 4, n), 4 * ones(1, degree)], degree);

% Optimization variables
t = sdpvar(1);
q = BSpline.sdpvar({B1, B2}, [1, 1]);
r = BSpline.sdpvar({B1, B2}, [2, 1]);

% parameter
th = parameter(2, 1);

% The optimization problem
p = -th(1) - th(2);
h1 = 2 * th(1) * th(1) * th(1) * th(1) - 8 * th(1) * th(1) * th(1) + 8 * th(1) * th(1) + 2 - th(2);
h2 = 4 * th(1) * th(1) * th(1) * th(1) - 32 * th(1) * th(1) * th(1) + 88 * th(1) * th(1) - 96 * th(1) + 36 - th(2);
sol = solvesdp([p - t - q - h1 * r(1) - h2 * r(2) == 0, q >= 0, r >= 0], -t, options);

mes = sprintf('Estimation: %0.5g \nLasserre2001: -5.5079', double(t));
disp(mes)

th1 = linspace(0, 3, 101);
th2 = linspace(0, 4, 101);
[T1, T2] = meshgrid(th1, th2);
figure
contour(T1, T2, -T1 - T2)
hold on
contour(T1, T2, 2 * T1.^4 - 8 * T1.^3 + 8 * T1.^2 + 2 - T2, [0, 0])
contour(T1, T2, 4 * T1.^4 - 32 * T1.^3 + 88 * T1.^2 - 96 * T1 + 36 - T2, [0, 0])
