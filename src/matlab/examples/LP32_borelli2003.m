clear all;
clc;

% Define Polygonal region
% =======================
vertices = [8 + 2 / 3, 8 + 1 / 3;
            1 + 2 / 3, -5 - 2 / 3;
            -5 - 1 / 3, 1 + 1 / 3;
            8 + 2 / 3, 8 + 1 / 3];

% Polygonal region
b_alpha = BSplineBasis([0, linspace(0, 1, 4), 1], 1);
c_alpha = mat2cell(vertices, [1, 1, 1, 1], 2);
alpha = BSpline(b_alpha, c_alpha);
x = linspace(0, 1, 100);
alpha_x = cell2mat(alpha.f(x));

% Add rays
b_beta = BSplineBasis([0, 0, 1, 1], 1);
c_beta = [0, 1];
beta = BSpline(b_beta, c_beta);

c_t1 = vertices(:, 1) * [0, 1];
c_t2 = vertices(:, 2) * [0, 1];
t1 = TensorSpline({b_alpha, b_beta}, c_t1);
t2 = TensorSpline({b_alpha, b_beta}, c_t2);

% Define optimization problem
degree = 1;
na = 11;
nb = 11;
Ba = BSplineBasis([0 * ones(1, degree) linspace(0, 1, na) ones(1, degree)], degree);
Bb = BSplineBasis([0 * ones(1, degree) linspace(0, 1, nb) ones(1, degree)], degree);

cvx_begin
    variable Ct1(length(Ba), length(Bb))
    variable Ct2(length(Ba), length(Bb))

    x = TensorSpline({Ba, Bb}, Ct1);
    y = TensorSpline({Ba, Bb}, Ct2);

    c1 = x + 3 * y + 2 * t1 - t2;
    c2 = 2 * x + y - t1 + 2 * t2;
    c3 = x - t1 - t2;
    objective = -2 * x - y;

    minimize(objective.integral)
    subject to
        c1.coeffs.coeffs2tensor <= 9
        c2.coeffs.coeffs2tensor <= 8
        c3.coeffs.coeffs2tensor <= 4
        Ct1 >= 0
        Ct2 >= 0
cvx_end

alpha = linspace(0, 1, 101);
beta = linspace(0, 1, 101);
T1 = t1.f({alpha, beta});
T2 = t2.f({alpha, beta});
x = TensorSpline({Ba, Bb}, Ct1);
y = TensorSpline({Ba, Bb}, Ct2);
objective = -2 * x - y;
surf(T1, T2, objective.f({alpha, beta}))
