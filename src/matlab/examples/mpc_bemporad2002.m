close all; clear all; clc;

% Problem data
% ============
% System
A = [0.7326, -0.0861; 0.1722 0.9909];
B = [0.0609; 0.0064];
C = [0, 1.4142];

% Control parameters
N_u = 2;  % Control horizon
N_y = 2;  % Prediction horizon
N_c = 1;  % Constraint horizon
% Q = eye(2);
% P = dlyap(A', Q);
% R = 0.01;
% bla
% A_pow_B = @(n) cell2mat(arrayfun(@(p) A^p * B, 0:n, 'uni', false));
% H = 2 * (A_pow_B(1)' * P * A_pow_B(1) + A_pow_B(0)' * Q * A_pow_B(0)) + R * eye(2);

% Hard coded, but should follow from P, Q, R
H = [1.5064, 0.4838; 0.4838, 1.5258];
F = [9.6652, 5.2115; 7.0732, -7.0879];
% H = [1.2152, 0.7848; 0.7848, 1.2152];
% F = [11.6581, 10.5542; 11.6925 10.5824];

% BSpline parameterization
degree = 2;
n_knots = 10;
knots = [-1.5 * ones(1, degree), linspace(-1.5, 1.5, n_knots), 1.5 * ones(1, degree)];
b = BSplineBasis(knots, degree);
x = Polynomial({[0; 0], [0; 1]; [1; 0], [0; 0]});
x = x.to_bspline({[-1.5, 1.5], [-1.5, 1.5]});

% Optimization problem
u = BSpline.sdpvar({b, b}, [N_u, 1]);
t = BSpline.sdpvar({b, b}, [1, 1]);
obj = 0.5 * t + x' * F * u;
LMI = [t, u'; u, inv(H)];  % Schur complement: makes it computationally cheaper for Yalmip
con = [u(1:N_c+1) <= 2, u(1:N_c+1) >= -2, LMI >= 0];
options = sdpsettings('verbose', 2, 'solver', 'sdpt3');
sol = solvesdp(con, obj.integral, options);

% Simulate MPC controller
u = double(u);
U = [];
xk = [1;1];
while norm(xk(:, end)) >= 1e-6
    U = [U, u(1).f(num2cell(xk(:, end)'))];
    xk = [xk, A * xk(:, end) + B * U(end)];
end
figure
plot(xk(1, :), xk(2, :))
xlabel('x_1')
ylabel('x_2')

t = 0.1 * (0:size(xk, 2)-1);
figure
subplot(211)
stairs(t, xk(1, :))
hold on
stairs(t, xk(2, :), 'r')
ylabel('states')
subplot(212)
stairs(t(1:end-1), U)
ylabel('input')
xlabel('t (s)')

% With additional state constraint
% ================================

% First solve phase 1
u = BSpline.sdpvar({b, b}, [N_u, 1]);
t = BSpline.sdpvar({b, b}, [1, 1]);
p = BSpline.sdpvar({b, b}, [1, 1]);
obj = p;
LMI = [t, u'; u, inv(H)];  % Schur complement: makes it computationally cheaper for Yalmip
con = [u(1:N_c+1) <= 2, u(1:N_c+1) >= -2, LMI >= 0, A * x + B * u(1) + p >= -0.5, A^2 * x + B * u(2) + A * B * u(1) + p >= -0.5];
options = sdpsettings('verbose', 1, 'solver', 'sdpt3');
sol = solvesdp(con, obj.integral, options);
p2 = double(p);
p_c = p2.coeffs.coeffs2tensor;
p_c(p_c <= 0) = 0;
p = BSpline({b, b}, p_c);

% Use solution of phase 1 to handle infeasibilities
u = BSpline.sdpvar({b, b}, [N_u, 1]);
t = BSpline.sdpvar({b, b}, [1, 1]);
obj = 0.5 * t + x' * F * u;
LMI = [t, u'; u, inv(H)];  % Schur complement: makes it computationally cheaper for Yalmip
con = [u(1:N_c+1) <= 2, u(1:N_c+1) >= -2, LMI >= 0, A * x + B * u(1) + p >= -0.5, A * x + B * u(2) + p >= -0.5];
options = sdpsettings('verbose', 1, 'solver', 'sdpt3');
sol = solvesdp(con, obj.integral, options);

% Simulate MPC controller
u = double(u);
U = [];
xk = [1;1];
while norm(xk(:, end)) >= 1e-6
    U = [U, u(1).f(num2cell(xk(:, end)'))];
    xk = [xk, A * xk(:, end) + B * U(end)];
end
figure
plot(xk(1, :), xk(2, :))
xlabel('x_1')
ylabel('x_2')
hold on
contour(-1:0.01:1, -1:0.01:1, p2.f({-1:0.01:1, -1:0.01:1})', [0, 0])


t = 0.1 * (0:size(xk, 2)-1);
figure
subplot(211)
stairs(t, xk(1, :))
hold on
stairs(t, xk(2, :), 'r')
ylabel('states')
subplot(212)
stairs(t(1:end-1), U)
ylabel('input')
xlabel('t (s)')