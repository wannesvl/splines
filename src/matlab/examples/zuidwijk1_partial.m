
% Test partial refinement of the approximations

clear all
close all
clear global
clc



% Inputs and Settings
% ===================

l = Polynomial([0, 1]);     % the parameter
L = 3;                      % range of l = [0, L]
d = 2;                      % degree
n = 11;                     % number of knots
% load zuidwijk1_partial_B2 
% B1 = B2; clear B2
B1 = BSplineBasis([0 * ones(1, d) , linspace(0, L, n) , L * ones(1, d)], d);
colors = colormap(lines(3));



% True solution
% =============

L_ = linspace(0, L, 1001)';
figure(1), hold on
plot(L_(L_<=1), 6 - 2 * L_(L_<=1), 'k')
plot(L_(L_<=2 & L_>=1), 4 ./ L_(L_<=2 & L_>=1), 'k')
plot(L_(L_>=2), 2 * ones(length(L_(L_ >= 2))), 'k')
xlabel('\theta'), ylabel('x_1(\theta) + 2 x_2(\theta)')



% Approximate solutions with basis B1
% ===================================

% Primal
x = BSpline.sdpvar(B1, [1, 2]);
obj = x(1) + 2 * x(2);
con = [x >= 0, x(2) <= 2, x(1) + l * x(2) <= 2];
solvesdp(con, -obj.integral, sdpsettings('verbose',1));
x1 = double(x);
f1 = double(obj);
plot(L_, f1.f(L_), 'Color', colors(1,:))

% Dual
y = BSpline.sdpvar(B1, [1, 2]);
obj = 2 * y(1) + 2 * y(2);
con = [y >= 0, y(2) >= 1, y(1) + l * y(2) >= 2];
solvesdp(con, obj.integral, sdpsettings('verbose',1));
y1 = double(y);
g1 = double(obj);
plot(L_, g1.f(L_), 'Color', colors(1,:))

% Duality gap
gap1 = g1 - f1;
Cgap1 = gap1.coeffs.coeffs2tensor;
B1_knots_ave = conv(ones(B1.degree,1) / B1.degree , B1.knots);
B1_knots_ave = B1_knots_ave(B1.degree + 1 : B1.degree + B1.length);
figure(2), hold all
plot(L_, gap1.f(L_), B1_knots_ave, Cgap1, ':o', 'Color', colors(1,:))
xlabel('\theta'), ylabel('duality gap')



% Knot insertion
% ==============

[~,imax] = max(Cgap1);                          % ??? which are the best knots to add ???
i = imax;
if abs(Cgap1(imax-1)-Cgap1(imax)) <= 0.15*Cgap1(imax)
    i = [imax-1; i];
end
if abs(Cgap1(imax+1)-Cgap1(imax)) <= 0.15*Cgap1(imax)
    i = [i; imax+1];
end
k2i = B1_knots_ave(i);
B2 = B1.insert_knots(k2i);
x1a = x1.insert_knots({k2i});
f1a = f1.insert_knots({k2i});
y1a = y1.insert_knots({k2i});
g1a = g1.insert_knots({k2i});

gap1a = g1a - f1a;
B2_knots_ave = conv(ones(B2.degree,1) / B2.degree , B2.knots);
B2_knots_ave = B2_knots_ave(B2.degree + 1 : B2.degree + B2.length);
figure(2), hold all
plot(B2_knots_ave,  gap1a.coeffs.coeffs2tensor, ':x', 'Color', colors(1,:))



% Full re-optimization after knot insertion
% =========================================

% Primal
x = BSpline.sdpvar(B2, [1, 2]);
obj = x(1) + 2 * x(2);
con = [x >= 0, x(2) <= 2, x(1) + l * x(2) <= 2];
solvesdp(con, -obj.integral, sdpsettings('verbose',1));
x2 = double(x);
f2 = double(obj);
figure(1), plot(L_, f2.f(L_), 'Color', colors(2,:))

% Dual
y = BSpline.sdpvar(B2, [1, 2]);
obj = 2 * y(1) + 2 * y(2);
con = [y >= 0, y(2) >= 1, y(1) + l * y(2) >= 2];
solvesdp(con, obj.integral, sdpsettings('verbose',1));
y2 = double(y);
g2 = double(obj);
figure(1), plot(L_, g2.f(L_), 'Color', colors(2,:))

% Duality gap
gap2 = g2 - f2;
figure(2)
plot(L_, gap2.f(L_), B2_knots_ave, gap2.coeffs.coeffs2tensor, ':x', 'Color', colors(2,:))



% Partial re-optimization after knot insertion
% ============================================

[~,i] = ismember(k2i, B2.knots);
i2opt = (min(i) - B2.degree - 1 : max(i) + B2.degree)';     % ??? which coeffs should be re-optimized ???

% Primal
Cx1a = x1a.coeffs.coeffs2tensor;
Cx2opt = sdpvar(length(i2opt), 2);
Cx = [Cx1a(1:i2opt(1)-1,:) ; Cx2opt; Cx1a(i2opt(end)+1:end,:)];
x = BSpline(B2, mat2cell(Cx, ones(B2.length,1),2));
obj = x(1) + 2 * x(2);
con = [x >= 0, x(2) <= 2, x(1) + l * x(2) <= 2];
solvesdp(con, -obj.integral, sdpsettings('verbose',1));
x3 = double(x);
f3 = double(obj);
figure(1), plot(L_, f3.f(L_), 'Color', colors(3,:))

% Dual 
Cy1a = y1a.coeffs.coeffs2tensor;
Cy2opt = sdpvar(length(i2opt), 2);
Cy = [Cy1a(1:i2opt(1)-1,:) ; Cy2opt; Cy1a(i2opt(end)+1:end,:)];
y = BSpline(B2, mat2cell(Cy, ones(B2.length,1),2));
obj = 2 * y(1) + 2 * y(2);
con = [y >= 0, y(2) >= 1, y(1) + l * y(2) >= 2];
solvesdp(con, obj.integral, sdpsettings('verbose',1));
y3 = double(y);
g3 = double(obj);
figure(1), plot(L_, g3.f(L_), 'Color', colors(3,:))

% Duality gap
gap3 = g3 - f3;
figure(2)
plot(L_, gap3.f(L_), B2_knots_ave, gap3.coeffs.coeffs2tensor, ':x', 'Color', colors(3,:))


save zuidwijk1_partial_B2 B2



