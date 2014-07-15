
clear all
close all
clear global
clc



% Inputs and settings
% ===================

% Problem data
A0  = [1 2 -3; 2 4 -1; -3 -1 3];
At1 = [1 -1 2; -1 1 3; 2 3 2];
At2 = [-1 1 0; 1 1 2; 0 2 -2];
Ax1 = [3 -2 4; -2 1 -2; 4 -2 -2];
Ax2 = [-3 1 1; 1 -2 -1; 1 -1 1];
Ax3 = [5 4 2; 4 1 1; 2 1 -1];
c = [1; -2; 1];

% Bspline parameterization
degree = 3;
n_knots = 21;
knots = [-2 * ones(1, degree), linspace(-2, 2, n_knots), 2 * ones(1, degree)];
Bth = BSplineBasis(knots, degree);
th = Polynomial({[0, 0], [0, 1]; [1, 0], [0, 0]});
th = th.to_bspline({[-2, 2], [-2, 2]});



% True solution by gridding
% =========================

try
    load SDP511_Bemporad2006
catch
    x = sdpvar(3, 1);
    t = sdpvar(2, 1);
    K = A0 + At1 * t(1) + At2 * t(2) + Ax1 * x(1) + Ax2 * x(2) + Ax3 * x(3);
    objective = x(1) - 2 * x(2) + x(3);
    prob = optimizer([K >= 0], objective, sdpsettings('verbose',0), t, x);

    th1 = linspace(-2, 2, 101)';
    th2 = linspace(-2, 2, 101)';
    [Th1, Th2] = meshgrid(th1, th2);
    F = zeros(size(Th1));
    for i1 = 1 : size(F,1)
        i1
        for i2 = 1 : size(F,2)
            [xk, errorcode] = prob{[Th1(i1,i2);Th2(i1,i2)]};
            F(i1,i2) = xk(1) - 2 * xk(2) + xk(3);
        end
    end
    
    save SDP511_Bemporad2006 th1 th2 Th1 Th2 F
end

F0 = F;
feas = zeros(size(F));
feas(isnan(F0)) = 1;
F(isnan(F0)) = 10;
surf(Th1, Th2, F), shading interp
xlabel('theta_1'), ylabel('theta_2'), zlabel('obj')



% Solve phase-1 problem (primal)
% ==============================

x = BSpline.sdpvar({Bth, Bth}, [3, 1]);
t = BSpline.sdpvar({Bth, Bth}, [1, 1]);
con = [A0 + At1 * th(1) + At2 * th(2) + Ax1 * x(1) + Ax2 * x(2) + Ax3 * x(3) + t * eye(3) >= 0 , ...
       t >= 0];
sol = solvesdp(con, t.integral, sdpsettings('verbose', 1, 'solver', 'sdpt3'));
t1 = double(t);

x = BSpline.sdpvar({Bth, Bth}, [3, 1]);
t = BSpline.sdpvar({Bth, Bth}, [1, 1]);
con = [A0 + At1 * th(1) + At2 * th(2) + Ax1 * x(1) + Ax2 * x(2) + Ax3 * x(3) + t * eye(3) >= 0];
sol = solvesdp(con, t.integral, sdpsettings('verbose', 1, 'solver', 'sdpt3'));
t2 = double(t);
t3 = BSpline({Bth, Bth}, max(t2.coeffs.coeffs2tensor, 0));

% Overview figure
T1 = t1.f({th1, th2})';
T2 = t2.f({th1, th2})';
T3 = t3.f({th1, th2})';
figure
subplot(231)
surf(Th1, Th2, T1), shading interp
xlabel('theta_1'), ylabel('theta_2'), zlabel('LMI relax')
subplot(232)
surf(Th1, Th2, T2), shading interp
xlabel('theta_1'), ylabel('theta_2'), zlabel('LMI relax')
subplot(233)
surf(Th1, Th2, T3), shading interp
xlabel('theta_1'), ylabel('theta_2'), zlabel('LMI relax')
subplot(234)
contour(Th1, Th2, feas, [0.5,0.5]), hold on
contour(Th1, Th2, T1, [1e-9,1e-6]), hold off
xlabel('theta_1'), ylabel('theta_2')
subplot(235)
contour(Th1, Th2, feas, [0.5,0.5]), hold on
contour(Th1, Th2, T2, [1e-9,1e-6]), hold off
xlabel('theta_1'), ylabel('theta_2')
subplot(236)
contour(Th1, Th2, feas, [0.5,0.5]), hold on
contour(Th1, Th2, T3, [1e-9,1e-6]), hold off
xlabel('theta_1'), ylabel('theta_2')



% Solve relaxed problems
% ======================

% Primal
t = t1;
x = BSpline.sdpvar({Bth, Bth}, [3, 1]);
con = [A0 + At1 * th(1) + At2 * th(2) + Ax1 * x(1) + Ax2 * x(2) + Ax3 * x(3) + t * eye(3) >= 0];
obj = c'*x;
sol = solvesdp(con, obj.integral); 
f1 = double(obj);

% Dual
Z = BSpline.sdpvar({Bth, Bth}, [3, 3], 'symmetric');
con = [trace(Z*Ax1) == c(1), trace(Z*Ax2) == c(2), trace(Z*Ax3) == c(3), Z >= 0];
obj = trace(Z * (A0 + At1*th(1) + At2*th(2) + t * eye(3)));
solvesdp(con, obj.integral, sdpsettings('verbose',1));
g1 = -double(obj);
g2 = -double(trace(Z * (A0 + At1*th(1) + At2*th(2))));

% Overview figure
F1 = f1.f({th1, th2})';     F1(isnan(F0)) = 10;
G1 = g1.f({th1, th2})';     G1(isnan(F0)) = 10;
G2 = g2.f({th1, th2})';     G2(isnan(F0)) = 10;
eF1 = abs(F1-F);            eF1u = max(eF1(:));
eG1 = abs(G1-F);            eG1u = max(eG1(:));
eG2 = abs(G2-F);            eG2u = max(eG2(:));

figure
subplot(231)
surf(Th1, Th2, F1), shading interp
xlabel('theta_1'), ylabel('theta_2'), zlabel('obj')
subplot(232)
surf(Th1, Th2, G1), shading interp
xlabel('theta_1'), ylabel('theta_2'), zlabel('obj')
subplot(233)
surf(Th1, Th2, G2), shading interp
xlabel('theta_1'), ylabel('theta_2'), zlabel('obj')
subplot(234)
surf(Th1, Th2, eF1, eF1 / max([eF1u,eG1u,eG2u])), shading interp
set(gca, 'CLim', [0,1], 'CLimMode', 'manual')
xlabel('theta_1'), ylabel('theta_2'), zlabel('e')
subplot(235)
surf(Th1, Th2, eG1, eG1 / max([eF1u,eG1u,eG2u])), shading interp
set(gca, 'CLim', [0,1], 'CLimMode', 'manual')
xlabel('theta_1'), ylabel('theta_2'), zlabel('e')
subplot(236)
surf(Th1, Th2, eG2, eG2 / max([eF1u,eG1u,eG2u])), shading interp
set(gca, 'CLim', [0,1], 'CLimMode', 'manual')
xlabel('theta_1'), ylabel('theta_2'), zlabel('e')




% TRASH
% % Alternatively relaxed dual
% Z = BSpline.sdpvar({Bth, Bth}, [3, 3], 'symmetric');
% con = [ trace(Z*Ax1) == c(1), trace(Z*Ax2) == c(2), trace(Z*Ax3) == c(3), ...
%         Z >= 0, ...
%         trace(Z * (A0 + At1*th(1) + At2*th(2))) >= -10 ];
% obj = trace(Z * (A0 + At1*th(1) + At2*th(2)));
% solvesdp(con, obj.integral, sdpsettings('verbose',1));
% g2 = -double(obj);
