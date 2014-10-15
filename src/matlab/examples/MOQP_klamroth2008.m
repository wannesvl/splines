
% Running example in [Klamroth, Miettinen, 2008], 
% adopted from [Schandl et al, 2001]

clear all 
close all
clear global
clc

addpath ../Function/
addpath ../Function/Coefficients/
addpath ../Function/Basis



% Problem Data
% ============

P1 = zeros(2);       q1 = [-1; -1];          r1 = 5;
P2 = 2/5 * eye(2);   q2 = 1/5 * [-10; -4];   r2 = 11/5;
A = [3, 1; 2, 1; 1, 2; -1, 0; 0, -1];
b = [-12; -9; -12; 0; 0];



% Feasible set and objective set
% ==============================

x1 = linspace(0,4,1e3)';
x2 = [12-3*x1, 9-2*x1, 6-0.5*x1];
xB = [ linspace(4,0,4e2)' , zeros(4e2,1)                  ;
       zeros(6e2,1)       , linspace(0,6,6e2)'            ;
       linspace(0,2,2e2)' , 6  - 0.5 * linspace(0,2,2e2)' ;
       linspace(2,3,1e2)' , 9  - 2   * linspace(2,3,1e2)' ;
       linspace(3,4,1e2)' , 12 - 3   * linspace(3,4,1e2)' ];
fB = [ 5 + xB*[-1; -1] , 1/5 * (11 + xB*[-10; -4] + (xB.^2)*[1; 1]) ];

figure(1)
area(xB(:,1), xB(:,2), 'FaceColor', [0.9,0.9,0.9], 'Edgecolor', 'none'), hold all
set(gca,'Layer','top')
plot(x1, x2, ':', 'Color', [0.7,0.7,0.7])
axis([0,4,0,6]), grid on
xlabel('x_1'), ylabel('x_2')

figure(2)
area(fB(:,1), fB(:,2), 'FaceColor', [0.9,0.9,0.9], 'Edgecolor', 'none'), hold all
set(gca,'Layer','top')
axis([-3,5,-4,5]), grid on
xlabel('f_1'), ylabel('f_2')



% Computing the trade-off curve by gridding
% =========================================

N = 100;
lam = linspace(0,1,N)';
try 
    load tradeoff_klamroth2008
catch
    x_opt = zeros(N,2);
    f_opt = zeros(N,2);
    g_opt = zeros(N,1);
    nu_opt = zeros(N,5);
    for i = 1:N
        x = sdpvar(2,1);
        constr = set( A*x + b <= 0 ); 
        f1 = 0.5 * x'*P1*x + q1'*x + r1;
        f2 = 0.5 * x'*P2*x + q2'*x + r2;
        goal = lam(i) * f1 + (1-lam(i)) * f2;
        sol = solvesdp(constr, goal, sdpsettings('verbose', 0));
        x_opt(i,:) = double(x)';
        f_opt(i,:) = [double(f1), double(f2)];

        P = lam(i) * P1 + (1-lam(i)) * P2;
        q = lam(i) * q1 + (1-lam(i)) * q2;
        r = lam(i) * r1 + (1-lam(i)) * r2;
        nu = sdpvar(5,1);
        t = sdpvar(1);
        constr = set([2*P, q+A'*nu; q'+nu'*A, t] >= 0) + set(nu >= 0);
        goal = t - r - nu'*b;
        sol = solvesdp(constr, goal, sdpsettings('verbose', 0));
        nu_opt(i,:) = double(nu)';
        g_opt(i) = -double(goal);
    end
    check_duality = max(abs( lam.*f_opt(:,1) + (1-lam).*f_opt(:,2) - g_opt ))
    save tradeoff_klamroth2008 x_opt f_opt nu_opt g_opt
end

figure(1), hold all
plot(x_opt(:,1), x_opt(:,2), 'k')
plot(x_opt(1,1), x_opt(1,2), 'ko')
plot(x_opt(end,1), x_opt(end,2), 'ko')

figure(2), hold all
plot(f_opt(:,1), f_opt(:,2), 'k')
plot(f_opt(1,1), f_opt(1,2), 'ko')
plot(f_opt(end,1), f_opt(end,2), 'ko')

figure(3)
plot(lam, [lam.*f_opt(:,1) + (1-lam).*f_opt(:,2), g_opt], 'k') 
xlabel('\lambda')
ylabel('f_{opt}')

figure(4)
plot(lam, x_opt(:,1), 'k'), hold all
plot(lam, x_opt(:,2), 'k:')
xlabel('\lambda')
ylabel('x_{opt}')



% Parameterized primal solution
% =============================

B1 = BSplineBasis([0 0 1 1], 1);
P = BSpline(B1, {P2; P1});
q = BSpline(B1, {q2; q1});
r = BSpline(B1, {r2; r1});

degree = 1;
Nknots = 5;
B = BSplineBasis([0 * ones(1, degree) linspace(0, 1, Nknots) ones(1, degree)], degree);

Cx = sdpvar(2, B.length, 'full');
x = BSpline(B, mat2cell(Cx, 2, ones(B.length,1)));
c = A*x + b;
f = 0.5 * x'*P*x + q'*x + r;
constr = set( c.coeffs.coeffs2tensor <= 0 );
goal = f.integral;
sol = solvesdp(constr, goal)
Cx = double(Cx);
goal_prim = double(goal)

x = BSpline(B, mat2cell(Cx, 2, ones(B.length,1)));
f = 0.5 * x'*P*x + q'*x + r;
f1 = 0.5 * x'*P1*x + q1'*x + r1;
f2 = 0.5 * x'*P2*x + q2'*x + r2;
lam = linspace(0,1,1e3)';
F = f.f(lam);
F1 = f1.f(lam);
F2 = f2.f(lam);
X = cell2mat(x.f(lam));
X1 = X(1:2:end);
X2 = X(2:2:end);

figure(1), hold all
plot(X1, X2, 'b')
% plot(Cx(1,:), Cx(2,:), 'rx')

figure(2), hold all
plot(F1, F2, 'b')
% plot(f1.coeffs.coeffs2tensor, f2.coeffs.coeffs2tensor, 'rx')

figure(3), hold all
plot(lam, F, 'b')

figure(4), hold all
plot(lam, X1, 'b')
plot(lam, X2, 'b:')


% % Reconstruct feasible dual solution (-> bad lower bound // is this ok to do??) 
% Cnu = mat2cell(dual(constr),5*ones(c.basis{1}.length,1),1);
% nu = BSpline(c.basis{1}, Cnu);
% Nu = nu.f({lam}).coeffs;
% g = zeros(size(lam));
% for i = 1:length(lam)
%     P = lam(i) * P1 + (1-lam(i)) * P2;
%     q = lam(i) * q1 + (1-lam(i)) * q2;
%     r = lam(i) * r1 + (1-lam(i)) * r2;
%     nu = Nu{i};
%     g(i) = r + nu' * b - 0.5 * (q+A'*nu)' * (P \ (q+A'*nu));
% end
% figure(3), hold all
% plot(lam, g, 'g')



% Parameterized dual solution
% ===========================

degree = 1;
Nknots = 5;
B = BSplineBasis([0 * ones(1, degree) linspace(0, 1, Nknots) ones(1, degree)], degree);

Cnu = sdpvar(5, B.length, 'full');
nu = BSpline(B, mat2cell(Cnu, 5, ones(B.length,1)));
Ct = sdpvar(1, B.length, 'full');
t = BSpline(B, Ct);
cLMI = [2*P, A'*nu+q; nu'*A+q', t];
cLMIc = cLMI.coeffs.coeffs2tensor;
constr = set(Cnu >= 0);
for i = 1:B.length
    constr = constr + set( cLMIc(3*(i-1)+1 : 3*i, :) >=0 );
end
g = t - r - nu'*b;
goal = g.integral;
sol = solvesdp(constr, goal)
goal_dual = double(goal)
Cnu = double(Cnu);      nu = BSpline(B, mat2cell(Cnu, 5, ones(B.length,1)));
Ct  = double(Ct);       t  = BSpline(B, Ct);
cLMI = [2*P, A'*nu+q; nu'*A+q', t];
g = t - r - nu'*b;
G = -g.f(lam);

figure(3), hold all
plot(lam, G, 'Color', [0,0.5,0])


% Get outer approximation of Pareto front
w = linspace(-2,2,1e2)';
figure(2), hold all
for i = 1:10:length(lam)
    c = [lam(i), 1-lam(i)];
    nc = [lam(i)-1, lam(i)];
    Fline = G(i)/norm(c)^2 * ones(size(w))*c + w*nc;
    plot(Fline(:,1), Fline(:,2), ':', 'Color', [0.5,0.8,0.5])
    plot(G(i)/norm(c)^2 * ones(size(w))*c(1), G(i)/norm(c)^2 * ones(size(w))*c(2), 'rx')
end


% Reconstruct feasible primal solution 
Bx = cLMI.basis{1};
intBx = Bx.integral;
Cx = zeros(3,Bx.length);
for i = 1:Bx.length-1
   V = dual(constr(i+1));
   V = V ./ intBx(i);
   [v,d] = eig(V);
   [~,im] = max(diag(d));
   Cx(:,i) = v(:,im) * sqrt(d(im,im));
end
i = Bx.length;
V = dual(constr(i+1));
V = V ./ intBx(i);
Cx(:,i) = V(:,end);

x = BSpline(Bx,  mat2cell(2*Cx(1:2,:), 2, ones(1,Bx.length)));
f = 0.5 * x'*P*x + q'*x + r;
f1 = 0.5 * x'*P1*x + q1'*x + r1;
f2 = 0.5 * x'*P2*x + q2'*x + r2;

F = f.f(lam);
F1 = f1.f(lam);
F2 = f2.f(lam);
X = cell2mat(x.f(lam));
X1 = X(1:2:end);
X2 = X(2:2:end);

figure(1), hold all
plot(X1, X2, 'Color', [0,0.5,0])

figure(2), hold all
plot(F1, F2, 'Color', [0,0.5,0])

figure(3), hold all
plot(lam, F, ':', 'Color', [0,0.5,0])

figure(4), hold all
plot(lam, X1, 'Color', [0,0.5,0])
plot(lam, X2, ':', 'Color', [0,0.5,0])
