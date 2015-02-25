close all; clear all; clc;

% ===========================================================================
% Example 1
% ===========================================================================

% Define Bspline basis
% ====================
degree = 3;
n_knots = 10;
Bth = BSplineBasis([-1.1, 1.1], degree, n_knots);
th = parameter(2);
t = BSpline.sdpvar({Bth, Bth}, [1, 1]);
K = [1 - 16 * th(1) * th(2), th(1); th(1), 1 - th(1) * th(1) - th(2) * th(2)];

% Solve phase-1 problem
% =====================
options = sdpsettings('verbose', 1, 'solver', 'sdpt3');
con = [K + t * speye(2) >= 0];
sol = solvesdp(con, t.integral, options);
t1 = value(t);

% Construct dual solution from lagrange multipliers
b_ = Bth.integral' * Bth.integral;
Y = {};
for i=1:length(con)
    Y{i} = dual(con(i)) / b_(i);
end
Y = reshape(Y, length(Bth), length(Bth));
Y = BSpline({Bth, Bth}, Y);
t1_ = trace(Y * K);

% Solve phase-1 problem with convexity constraint
% ===============================================
g = t.gradient;
con = [K + t * speye(2) >= 0, t.hessian >= 0];
sol = solvesdp(con, t.integral, options);
t2 = value(t);

% Construct dual solution from lagrange multipliers
Y = {};
for i=1:length(Bth)^2
    Y{i} = dual(con(i)) / b_(i);
end
Y = reshape(Y, length(Bth), length(Bth));
Y = BSpline({Bth, Bth}, Y);
t2_ = trace(Y * K);

% Solve phase-1 dual
% ==================
Z = BSpline.sdpvar({Bth, Bth}, [2, 2], 'symmetric');
obj = trace(Z * K);

con = [Z >= 0, trace(Z) == 1];
sol = solvesdp(con, obj.integral, options);
t3 = -value(obj);
 
% Solve phase-1 dual convex
% =========================
sol = solvesdp([Z >= 0, trace(Z) == 1, obj.hessian <= 0], obj.integral, options);
t4 = -value(obj);
 
% Plot inner and outer approximations
% ===================================
th1 = linspace(-1.1, 1.1, 201);
th2 = linspace(-1.1, 1.1, 201);
[T1, T2] = meshgrid(th1, th2);
F = 1 - 2*T1.^2 - T2.^2 - 16 * (T1 .* T2 - T1.^3 .* T2 - T1 .* T2 .^3);
F(F > 0 & 1 - T1.^2 - T2.^2 < 0) = -20;

figure
subplot(121)
contour(T1, T2, t1.f({th1, th2})', [0, 0])
hold all
contour(T1, T2, t3.f({th1, th2})', [0, 0])
contour(T1, T2, t1_.f({th1, th2})', [0, 0], 'k')
contour(T1, T2, F, [0, 0], 'r')
axis equal
subplot(122)
hold all
contour(T1, T2, t2.f({th1, th2})', [0, 0])
contour(T1, T2, t4.f({th1, th2})', [0, 0])
contour(T1, T2, t2_.f({th1, th2})', [0, 0], 'k')
contour(T1, T2, F, [0, 0], 'r')
axis equal

% ===========================================================================
% Example 2
% ===========================================================================
n_knots = 8;  % Reduce knots to limit computation times
Bth = BSplineBasis([-2, 2], degree, n_knots);
th = parameter(3);
t = BSpline.sdpvar({Bth, Bth, Bth}, [1, 1]);
K = [1 - th(3) * th(3), th(1) - th(2) * th(3), th(2) - th(1) * th(3); 
     th(1) - th(2) * th(3), 1 + th(1) * th(1) - th(2) * th(2) - th(3) * th(3), th(1) - th(2) * th(3);
     th(2) - th(1) * th(3), th(1) - th(2) * th(3), 1 - th(3) * th(3)];
options = sdpsettings('verbose', 1, 'solver', 'sdpt3');

% Solve phase-1 problem
% =====================
sol = solvesdp([K + t * speye(3) >= 0], t.integral, options);
t1 = value(t);

figure
hold on
th1 = -2:0.1:2;
[T1, T2, T3] = meshgrid(th1, th1, th1);
v = permute(t1.f({th1, th1, th1}), [2, 1, 3]);
p1 = patch([-3, -1, 1, -3], [3, -1, -1, 3], [-1, 1, -1, -1], 'b');
p2 = patch([-1, 1, 3, -1], [-1, -1, 3, -1], [1, -1, 1, 1], 'b');
[th1, th3] = meshgrid(-3:0.05:1, -1:0.05:1);
th1 = bsxfun(@plus, th1, (-1:0.05:1)' + 1);
p3 = surf(th1, -th3.^2 + th1 .* th3 + 1, th3);
set(p3,'FaceColor','blue','EdgeColor','none');
p = patch(isosurface(T1, T2, T3, v, 0));
isonormals(T1, T2, T3, v, p);
set(p,'FaceColor','red','EdgeColor','none');
daspect([1, 1, 1])
xlabel('th1'); ylabel('th2'); zlabel('th3');
alpha(p1, 0.2);
alpha(p2, 0.2);
alpha(p3, 0.2);
view(3); axis([-3 3 -1 3 -1 1 0 1]);
camlight 
lighting gouraud

% Solve phase-1 problem with convexity constraint
% ===============================================
sol = solvesdp([K + t * eye(3) >= 0, t.hessian >= 0], t.integral, options);
t2 = value(t);

figure
hold on
th1 = -2:0.05:2;
[T1, T2, T3] = meshgrid(th1, th1, th1);
v = permute(t2.f({th1, th1, th1}), [2, 1, 3]);
p1 = patch([-3, -1, 1, -3], [3, -1, -1, 3], [-1, 1, -1, -1], 'b');
p2 = patch([-1, 1, 3, -1], [-1, -1, 3, -1], [1, -1, 1, 1], [0,0,1,0.2]);
[th1, th3] = meshgrid(-3:0.05:1, -1:0.05:1);
th1 = bsxfun(@plus, th1, (-1:0.05:1)' + 1);
p3 = surf(th1, -th3.^2 + th1 .* th3 + 1, th3);
set(p3,'FaceColor','blue','EdgeColor','none');
p = patch(isosurface(T1, T2, T3, v, 0));
isonormals(T1, T2, T3, v, p);
set(p,'FaceColor','red','EdgeColor','none');
daspect([1, 1, 1])
xlabel('th1'); ylabel('th2'); zlabel('th3');
alpha(p1, 0.2);
alpha(p2, 0.2);
alpha(p3, 0.2);
view(3); axis([-3 3 -1 3 -1 1 0 1]);
camlight 
lighting gouraud

% Solve phase-1 dual
% ==================
Bth = BSplineBasis([-3, 3], degree, n_knots);
Z = BSpline.sdpvar({Bth, Bth, Bth}, [3, 3], 'symmetric');
obj = trace(Z * K);
sol = solvesdp([Z >= 0, trace(Z) == 1], obj.integral, options);
t3 = -value(obj);

figure
hold on
th1 = -3:0.05:3;
[T1, T2, T3] = meshgrid(th1, th1, th1);
v = permute(t3.f({th1, th1, th1}), [2, 1, 3]);
p1 = patch([-3, -1, 1, -3], [3, -1, -1, 3], [-1, 1, -1, -1], 'r');
p2 = patch([-1, 1, 3, -1], [-1, -1, 3, -1], [1, -1, 1, 1], 'r');
[th1, th3] = meshgrid(-3:0.05:1, -1:0.05:1);
th1 = bsxfun(@plus, th1, (-1:0.05:1)' + 1);
p3 = surf(th1, -th3.^2 + th1 .* th3 + 1, th3);
set(p3,'FaceColor','red','EdgeColor','none');
p = patch(isosurface(T1, T2, T3, v, 0));
isonormals(T1, T2, T3, v, p);
set(p,'FaceColor','blue','EdgeColor','none');
daspect([1, 1, 1])
xlabel('th1'); ylabel('th2'); zlabel('th3');
alpha(p, 0.5);
alpha(p1, 0.2);
alpha(p2, 0.2);
alpha(p3, 0.2);
view(3); axis([-3 3 -1 3 -1 1 0 1]);
camlight 
lighting gouraud
