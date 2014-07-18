
% Test different norms

clear all
close all
clear global
clc



% Inputs and Settings
% ===================

l = Polynomial([0, 1]);     % the parameter
L = 3;                      % range of l = [0, L]
d = 2;                      % degree
n = 21;                     % number of knots
B = BSplineBasis([0 * ones(1, d) , linspace(0, L, n) , L * ones(1, d)], d);
colors = colormap(lines(3));



% True solution
% =============

L_ = linspace(0, L, 1001)';
figure(1), hold on
plot(L_(L_<=1), 6 - 2 * L_(L_<=1), 'k')
plot(L_(L_<=2 & L_>=1), 4 ./ L_(L_<=2 & L_>=1), 'k')
plot(L_(L_>=2), 2 * ones(length(L_(L_ >= 2))), 'k')
xlabel('\theta'), ylabel('x_1(\theta) + 2 x_2(\theta)')



% Approximate solutions with L1 norm
% ==================================

% Primal
x = BSpline.sdpvar(B, [1, 2]);
obj = x(1) + 2 * x(2);
con = [x >= 0, x(2) <= 2, x(1) + l * x(2) <= 2];
solvesdp(con, -obj.integral, sdpsettings('verbose',1));
x1 = double(x);
f1 = double(obj);
plot(L_, f1.f(L_), 'Color', colors(1,:))

% Dual
y = BSpline.sdpvar(B, [1, 2]);
obj = 2 * y(1) + 2 * y(2);
con = [y >= 0, y(2) >= 1, y(1) + l * y(2) >= 2];
solvesdp(con, obj.integral, sdpsettings('verbose',1));
y1 = double(y);
g1 = double(obj);
plot(L_, g1.f(L_), 'Color', colors(1,:))

% Duality gap
gap1 = g1 - f1;
B_knots_ave = conv(ones(B.degree,1) / B.degree , B.knots);
B_knots_ave = B_knots_ave(B.degree + 1 : B.degree + B.length);
figure(2), hold all
plot(L_, gap1.f(L_), B_knots_ave, gap1.coeffs.coeffs2tensor, ':x', 'Color', colors(1,:))
xlabel('\theta'), ylabel('duality gap')

C1 = [x1(2)-2; x1(1)+l*x1(2)-2; -x1'];
Y1 = [y1'; -1+y1(2); -2+y1(1)+l*y1(2)];
s = -C1'*Y1;
plot(L_, s.f(L_), 'k')



% Approximate solutions with L2 norm
% ==================================

% Primal + dual
x = BSpline.sdpvar(B, [1, 2]);
y = BSpline.sdpvar(B, [1, 2]);
con = [ x >= 0, x(2) <= 2, x(1) + l * x(2) <= 2, ...
        y >= 0, y(2) >= 1, y(1) + l * y(2) >= 2 ];
obj = ( (2 * y(1) + 2 * y(2)) - (x(1) + 2 * x(2)) ) * ( (2 * y(1) + 2 * y(2)) - (x(1) + 2 * x(2)) );
solvesdp(con, obj.integral, sdpsettings('verbose',1));
x2 = double(x);
f2 = x2(1) + 2 * x2(2);
y2 = double(y);
g2 = 2 * y2(1) + 2 * y2(2);
figure(1), plot(L_, f2.f(L_), 'Color', colors(2,:))
figure(1), plot(L_, g2.f(L_), 'Color', colors(2,:))

% Duality gap
gap2 = g2 - f2;
figure(2)
plot(L_, gap2.f(L_), B_knots_ave, gap2.coeffs.coeffs2tensor, ':x', 'Color', colors(2,:))



% Approximate solutions with Linfinity norm
% =========================================

% Primal + dual
x = BSpline.sdpvar(B, [1, 2]);
y = BSpline.sdpvar(B, [1, 2]);
t = sdpvar(1);
con = [ x >= 0, x(2) <= 2, x(1) + l * x(2) <= 2, ...
        y >= 0, y(2) >= 1, y(1) + l * y(2) >= 2, ...
        (2 * y(1) + 2 * y(2)) - (x(1) + 2 * x(2)) <= t ];
solvesdp(con, t, sdpsettings('verbose',1));
x3 = double(x);
f3 = x3(1) + 2 * x3(2);
y3 = double(y);
g3 = 2 * y3(1) + 2 * y3(2);
figure(1), plot(L_, f3.f(L_), 'Color', colors(3,:))
figure(1), plot(L_, g3.f(L_), 'Color', colors(3,:))

% Duality gap
gap3 = g3 - f3;
figure(2)
plot(L_, gap3.f(L_), B_knots_ave, gap3.coeffs.coeffs2tensor, ':x', 'Color', colors(3,:))
