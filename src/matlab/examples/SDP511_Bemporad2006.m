close all; clear all; clc;

% Problem data
% ============
A0 = [1 2 -3; 2 4 -1; -3 -1 3];
At1 = [1 -1 2; -1 1 3; 2 3 2];
At2 = [-1 1 0; 1 1 2; 0 2 -2];
Ax1 = [3 -2 4; -2 1 -2; 4 -2 -2];
Ax2 = [-3 1 1; 1 -2 -1; 1 -1 1];
Ax3 = [5 4 2; 4 1 1; 2 1 -1];


% Feasible region
% ===============
degree = 3;
knots = [0 * ones(1, degree), linspace(0, 2 * pi, 20), 2 * pi * ones(1, degree)];
step = 1 / 50;
knots = 2 * pi * ((-degree * step):step:(1 + degree * step));
b_phi = BSplineBasis(knots, degree);

% Define objective integral
[b1, P1] = b_phi.derivative(1);
[b2, P2] = b_phi.derivative(2);
[b3, P3] = b_phi.derivative(3);
f_sin = @(x, c) -b_phi.f(x) * c .* cos(x) + b1.f(x) * P1 * c .* sin(x) + b2.f(x) * P2 * c .* cos(x)- b3.f(x) * P3 * c .* sin(x);
f_cos = @(x, c) b_phi.f(x) * c .* sin(x) + b1.f(x) * P1 * c .* cos(x) - b2.f(x) * P2 * c .* sin(x) - b3.f(x) * P3 * c .* cos(x);
uk = unique(knots)';
eps = 1e-15;
integ = @(f, c) sum(f(uk(2:end) - eps, c) - f(uk(1:end - 1) + eps, c));

cvx_solver sdpt3
cvx_begin sdp
    variable Ct1(length(b_phi))
    variable Ct2(length(b_phi))
    variable Cx1(length(b_phi))
    variable Cx2(length(b_phi))
    variable Cx3(length(b_phi))

    x1 = BSpline(b_phi, Cx1);
    x2 = BSpline(b_phi, Cx2);
    x3 = BSpline(b_phi, Cx3);
    t1 = BSpline(b_phi, Ct1);
    t2 = BSpline(b_phi, Ct2);

    K = A0 + At1 * t1 + At2 * t2 + Ax1 * x1 + Ax2 * x2 + Ax3 * x3;

    maximize(integ(f_sin, Ct1) + integ(f_cos, Ct2))  % Add regularization for spacing points evenly
    subject to
        for i=1:length(K.coeffs)
            K.coeffs.coeffs{i} >= 0;
        end
        Ct2 >= -2;
        Ct1 <= 2;
        % Periodicity constraints
        for i=1:degree
            Ct1(i) == Ct1(end - degree + i);
            Ct2(i) == Ct2(end - degree + i);
        end
        % for i=1:1
        %     Ct1(i) == Ct1(end - i + 1);
        %     Ct2(i) == Ct2(end - i + 1);
        % end
cvx_end

t1 = BSpline(b_phi, Ct1);
t2 = BSpline(b_phi, Ct2);

figure
t = linspace(0, 2 * pi, 501);
plot(t1.f(t), t2.f(t))

% Solve original problem
% ======================

% fill boundary
b_beta = BSplineBasis([0, 0, 1, 1], 1);
c_beta = [0, 1];
beta = BSpline(b_beta, c_beta);

c_t1 = (Ct1 - 2) * c_beta + 2;
c_t2 = (Ct2 + 2) * c_beta - 2;
t1 = BSpline({b_phi, b_beta}, c_t1);
t2 = BSpline({b_phi, b_beta}, c_t2);

% Define optimization problem
degree = 3;
na = 5;
nb = 20;
Ba = b_phi;
Bb = BSplineBasis([0 * ones(1, degree), linspace(0, 1, nb), ones(1, degree)], degree);

cvx_begin sdp
    variable Cx1(length(Ba), length(Bb))
    variable Cx2(length(Ba), length(Bb))
    variable Cx3(length(Ba), length(Bb))

    x1 = BSpline({Ba, Bb}, Cx1);
    x2 = BSpline({Ba, Bb}, Cx2);
    x3 = BSpline({Ba, Bb}, Cx3);

    K = A0 + At1 * t1 + At2 * t2 + Ax1 * x1 + Ax2 * x2 + Ax3 * x3;
    objective = x1 - 2 * x2 + x3;

    minimize(objective.integral)
    subject to
        for i=1:numel(K.coeffs.coeffs)
            K.coeffs.coeffs{i} >= 0;
        end
cvx_end

a = linspace(0, 2 * pi, 101);
b = linspace(0, 1, 101);
T1 = t1.f({a, b});
T2 = t2.f({a, b});
x1 = BSpline({Ba, Bb}, Cx1);
x2 = BSpline({Ba, Bb}, Cx2);
x3 = BSpline({Ba, Bb}, Cx3);
objective = x1 - 2 * x2 + x3;
O_approx = objective.f({a, b});
figure
surf(T1, T2, O_approx, 'EdgeColor', 'none')
camlight left; light; lighting phong; shading interp; alpha(0.5);

% Solve gridded problem
% =====================

O = zeros(21, 21);
for i=1:5:length(a)
    for j=1:5:length(b)
        cvx_begin sdp quiet
            variable x1
            variable x2
            variable x3

            K = A0 + At1 * T1(i, j) + At2 * T2(i, j) + Ax1 * x1 + Ax2 * x2 + Ax3 * x3;
            objective = x1 - 2 * x2 + x3;

            minimize(objective)
            subject to
                K >= 0;
        cvx_end
        O(ceil(i / 5), ceil(j / 5)) = x1 - 2 * x2 + x3;
    end
end

hold on;
surf(T1(1:5:length(a), 1:5:length(b)), T2(1:5:length(a), 1:5:length(b)), O)
camlight left; light; lighting phong; shading interp; alpha(0.5);


figure
surf(T1(1:5:length(a), 1:5:length(b)), T2(1:5:length(a), 1:5:length(b)), O - O_approx(1:5:length(a), 1:5:length(b)))