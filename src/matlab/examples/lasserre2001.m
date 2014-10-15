close all; clear all; clc;
options = sdpsettings('verbose', 1, 'solver', 'gurobi');
% Polynomial optimization Example Lasserre2001

% bases
degree = 2;
n = 5;
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
sol = solvesdp([p - t - h1 * r(1) - h2 * r(2) >= 0, r >= 0], -t);


x = fsolve(@(x) 2 * x.^4 - 8 * x.^3 + 8 * x.^2 + 2 -(4 * x.^4 - 32 * x.^3 + 88 * x.^2 - 96 * x + 36), 2.3, optimoptions('fsolve', 'TolX', 1e-10, 'TolFun', 1e-10));
y = 2 * x.^4 - 8 * x.^3 + 8 * x.^2 + 2;

mes = sprintf('Estimation: %0.5g \nLasserre2001: -5.5079\nActual optimum: %0.5g', double(t), -x - y);
disp(mes)

th1 = linspace(0, 3, 501);
th2 = linspace(0, 4, 501);
[T1, T2] = meshgrid(th1, th2);
Z = -T1 - T2;
figure
[ch1, ~] = contour(T1, T2, 2 * T1.^4 - 8 * T1.^3 + 8 * T1.^2 + 2 - T2, [0, 0]);
[ch2, ~] = contour(T1, T2, 4 * T1.^4 - 32 * T1.^3 + 88 * T1.^2 - 96 * T1 + 36 - T2, [0, 0]);
clf
ch1 = sortrows(ch1(:, 2:end)', 1);
ch2 = sortrows(ch2', 1);
yh1 = interp1(ch1(:, 1), ch1(:, 2), th1, 'linear', nan);
yh2 = interp1(ch2(:, 1), ch2(:, 2), th1, 'linear', nan);
plot(th1, min(yh1, yh2))
Z(T2 >= repmat(min(yh1, yh2), 501, 1)) = nan;
hold on
contourf(T1, T2, Z);
xlabel('x_1')
ylabel('x_2')
axis([0 3 0 4])
