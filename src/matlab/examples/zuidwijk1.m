% close all
clear all
clc

d = 3;
L = 4;
n = 21;
Bl = BSplineBasis([0 * ones(1, d) linspace(0, L, n) L * ones(1, d)], d);
l = BSpline(BSplineBasis([0, 0, L, L], 1), [0, L]');

cvx_solver sedumi  % sdpt3 fails on dual for large n
% Primal problem
cvx_begin
    variable c1(length(Bl))
    variable c2(length(Bl))

    x1 = BSpline(Bl, c1);
    x2 = BSpline(Bl, c2);
    obj = x1 + 2 * x2;
    con = [x1, x2, x1 + l * x2];

    minimize (-obj.integral)
    subject to
        con(1).coeffs >= 0;
        con(2).coeffs >= 0;
        con(2).coeffs <= 2;
        con(3).coeffs <= 2;
cvx_end

x = linspace(0, L, 101);
x1 = BSpline(Bl, c1);
x2 = BSpline(Bl, c2);
obj = x1 + 2 * x2;
plot(x, obj.f(x))
hold on
plot(x(x<=1), 6 - 2 * x(x<=1), 'r')
plot(x(x<=2 & x>=1), 4 ./ x(x<=2 & x>=1), 'r')
plot(x(x>=2), 2 * ones(length(x(x >= 2))), 'r')

% Dual problem
cvx_begin
    variable c1(length(Bl))
    variable c2(length(Bl))
    variable c3(length(Bl))
    variable c4(length(Bl))

    y1 = BSpline(Bl, c1);
    y2 = BSpline(Bl, c2);
    y3 = BSpline(Bl, c3);
    y4 = BSpline(Bl, c4);
    obj = 2 * y3 + 2 * y4;
    con = [y4 - y1, -y2 + y3 + l * y4];

    minimize (obj.integral)
    subject to
        con(1).coeffs == 1;
        con(2).coeffs == 2;
        y1.coeffs >= 0;
        y2.coeffs >= 0;
        y3.coeffs >= 0;
        y4.coeffs >= 0;
cvx_end
y3 = BSpline(Bl, c3);
y4 = BSpline(Bl, c4);
obj = 2 * y3 + 2 * y4;
plot(x, obj.f(x), 'g')

