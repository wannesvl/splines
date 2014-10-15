close all; clear all; clc;

% Problem data
% ============
A0 = [1 2 -3; 2 4 -1; -3 -1 3];
At1 = [1 -1 2; -1 1 3; 2 3 2];
At2 = [-1 1 0; 1 1 2; 0 2 -2];
Ax1 = [3 -2 4; -2 1 -2; 4 -2 -2];
Ax2 = [-3 1 1; 1 -2 -1; 1 -1 1];
Ax3 = [5 4 2; 4 1 1; 2 1 -1];

% Define Bspline basis
% ====================
degree = 3;
n_knots = 20;
knots = [-2 * ones(1, degree), linspace(-2, 2, n_knots), 2 * ones(1, degree)];
Bth = BSplineBasis(knots, degree);
th = Polynomial({[0, 0], [0, 1]; [1, 0], [0, 0]});

% Solve phase-1 problem
% =====================
x = BSpline.sdpvar({Bth, Bth}, [1, 3]);
t = BSpline.sdpvar({Bth, Bth}, [1, 1]);
obj = x(1) - 2 * x(2) + x(3);
K = A0 + At1 * th(1) + At2 * th(2) + Ax1 * x(1) + Ax2 * x(2) + Ax3 * x(3);
options = sdpsettings('verbose', 1, 'solver', 'sdpt3');
sol = solvesdp([K + t * eye(3) >= 0, t >= 0], t.integral, options);  % Phase 1
t = double(t);
sol = solvesdp([K + t * eye(3) >= 0], obj.integral, options);  % True problem

% Plot results
% ============
t = double(t);
obj = double(obj);
th1 = linspace(-2, 2, 101);
th2 = linspace(-2, 2, 101);
[T1, T2] = meshgrid(th1, th2);
figure
surf(T1, T2, t.f({th1, th2})')
camlight left; light; lighting phong; shading interp;

figure
surf(T1, T2, obj.f({th1, th2})')
camlight left; light; lighting phong; shading interp; alpha(0.5);

% Solve gridded problem
% =====================
try
    load SDP511_Bemporad2006_phase1
catch
    O = zeros(21, 21);
    x = sdpvar(3, 1);
    t = sdpvar(2, 1);
    K = A0 + At1 * t(1) + At2 * t(2) + Ax1 * x(1) + Ax2 * x(2) + Ax3 * x(3);
    objective = x(1) - 2 * x(2) + x(3);
    prob = optimizer([K >= 0], objective, options, t, x);
    for i=1:5:length(th1)
        for j=1:5:length(th2)
            [xk, errorcode] = prob{[th1(i); th2(j)]};
            O(ceil(i / 5), ceil(j / 5)) = xk(1) - 2 * xk(2) + xk(3);
        end
    end
    T1_ = T1(1:5:end, 1:5:end);
    T2_ = T2(1:5:end, 1:5:end);
    save('SDP511_Bemporad2006_phase1', 'O', 'T1_', 'T2_');
end

hold on;
surf(T1_, T2_, O')
camlight left; light; lighting phong; shading interp; alpha(0.5);

% Maximum error
err= obj.f({th1(1:5:end), th2(1:5:end)}) - O;
max(err(:))
figure
surf(T1_, T2_, err')
camlight left; light; lighting phong; shading interp;
