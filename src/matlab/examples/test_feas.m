
% Getting a spline bounded approximation of an ellipsoid

clear all
close all
clear global
clc

addpath ../Function/
addpath ../Function/Coefficients/
addpath ../Function/Basis



% The ellipsoid
% =============

a = pi/4;
U = [ cos(a), -sin(a);
      sin(a),  cos(a)];
D = diag([4;2]);
P = U*D*U';     
P2 = P*P;
Pinv = inv(P);
phi = linspace(0,2*pi,1e5)';
x = P(1,:) * [cos(phi'); sin(phi')];  x = x(:);
y = P(2,:) * [cos(phi'); sin(phi')];  y = y(:);

figure(1)
plot(x, y, 'k')
xlabel('x')
ylabel('y')



% Preprocessing : Get B-spline approx of cos(phi) and sin(phi)
% =============

Np = 1e2;
p = linspace(0,1,Np)';
phi = 2*pi*p;
degree = 2;
Nknots = 8; 
b = BSplineBasis([0 * ones(1, degree) linspace(0, 1, Nknots) ones(1, degree)] , degree); 
B = b.f(p);
c = B \ [cos(phi), sin(phi)];
bcos = BSpline(b, c(:,1));
bsin = BSpline(b, c(:,2));
bdir = [bcos; bsin];
bNdir = [-bsin; bcos];

% figure(100)
% plot(phi, [cos(phi), bcos.f(p), sin(phi), bsin.f(p)])
% figure(101)
% plot(cos(phi), sin(phi), bcos.f(p), bsin.f(p))



% True parameterizations in terms of p
% ====================================

try 
    eval(['load pars_ellipsoid_Np', num2str(Np)])
    
catch
    
    % Directional minimization
    f1 = zeros(Np,1);
    X1 = zeros(Np,2);
    for i = 1:Np
        x = sdpvar(2,1);
        constr = set([P2, x; x', 1] >= 0);
        goal = -[cos(phi(i)), sin(phi(i))] * x;
        sol = solvesdp(constr, goal, sdpsettings('solver', 'sdpt3', 'verbose', 0));
        X1(i,:) = double(x)';
        f1(i) = -double(goal);
    end

    % Radial expansion
    f2 = zeros(Np,1);
    X2 = zeros(Np,2);
    for i = 1:Np
        d = [cos(phi(i)); sin(phi(i))];
        r = sdpvar(1,1);
        constr = set([P2, r*d; r*d', 1] >= 0);
        goal = -r;
        sol = solvesdp(constr, goal, sdpsettings('solver', 'sdpt3', 'verbose', 0));
        X2(i,:) = double(r)*d';
        f2(i) = -double(goal);
    end

    % Saving the results
    eval(['save pars_ellipsoid_Np', num2str(Np), ' f1 X1 f2 X2'])
end


figure(1), hold all
plot(X1(:,1), X1(:,2), '.')
plot(X2(:,1), X2(:,2), '.')

figure(2)
plot(p, f1, p, f2)
xlabel('\phi / 2\pi')
ylabel('f')



% Approximate primals
% ===================

degree = 1;
Nknots = 11;
B = BSplineBasis([0 * ones(1, degree), linspace(0, 1, Nknots), 1 * ones(1, degree)]' , degree); 


% Directional minimization
Cx = sdpvar(2, B.length, 'full');
x = BSpline(B, mat2cell(Cx, 2, ones(B.length,1)));
constr = set([P2, x; x', 1] >= 0);   
for i = 1:degree
    dx = x.derivative(i-1, 1);
    constr = constr + set(dx.coeffs.coeffs{1} == dx.coeffs.coeffs{end});
end
f = - bdir' * x;
sol = solvesdp(constr, f.integral, sdpsettings('solver', 'sdpt3'))
Cx1h = double(Cx);          x1h = BSpline(B, mat2cell(Cx1h, 2, ones(B.length,1)));
X1h = cell2mat(x1h.f(p));   X1h = [X1h(1:2:end), X1h(2:2:end)];
f1h = bdir' * x1h;

figure(1), hold all
plot(X1h(:,1), X1h(:,2), '.-')

figure(2), hold all
plot(p, f1h.f(p))


% % Radial expansion    <??? first try, yet not a very good approach ???>
% Cr = sdpvar(1, B.length, 'full');
% r = BSpline(B, Cr);
% constr = constr + set([P2, r*bdir; r*bdir', 1] >= 0);   
% for i = 1:degree
%     dr = r.derivative(i-1, 1);
%     constr = constr + set(dr.coeffs.coeffs{1} == dr.coeffs.coeffs{end});
% end
% sol = solvesdp(constr, -r.integral, sdpsettings('solver', 'sdpt3'))
% Cr = double(Cr);
% r = BSpline(B, Cr);
% x2h = r*bdir;
% X2h = cell2mat(x2h.f(p));
% X2h = [X2h(1:2:end), X2h(2:2:end)];
% 
% figure(1), hold all
% plot(X2h(:,1), X2h(:,2), '.-')
% 
% figure(2), hold all
% plot(p, r.f(p))


% Radial expansion
p_knots = linspace(0, 1, Nknots);
Cx = sdpvar(2, B.length, 'full');
x = BSpline(B, mat2cell(Cx, 2, ones(B.length,1)));
constr = set([P2, x; x', 1] >= 0);
xf = x.f(p_knots);
Nf = bNdir.f(p_knots);
for i = 1:length(p_knots)-1
    constr = constr + set(xf{i}' * Nf{i} == 0);   
end
for i = 1:degree
    dx = x.derivative(i-1, 1);
    constr = constr + set(dx.coeffs.coeffs{1} == dx.coeffs.coeffs{end});
end
f = - bdir' * x;
sol = solvesdp(constr, f.integral, sdpsettings('solver', 'sdpt3'))
Cx2h = double(Cx);          x2h = BSpline(B, mat2cell(Cx2h, 2, ones(B.length,1)));
X2h = cell2mat(x2h.f(p));   X2h = [X2h(1:2:end), X2h(2:2:end)];
f2h = bdir' * x2h;

figure(1), hold all
plot(X2h(:,1), X2h(:,2), '.-')

figure(2), hold all
plot(p, f2h.f(p))



% Approximate duals
% =================

degree = 2;
Nknots = 10;
B = BSplineBasis([0 * ones(1, degree), linspace(0, 1, Nknots), 1 * ones(1, degree)]' , degree); 


% Dual for directional minimization     <??? How to use splinevar here ???>
cV11 = sdpvar(B.length,1);  V11 = BSpline(B, cV11);
cV12 = sdpvar(B.length,1);  V12 = BSpline(B, cV12);
cV22 = sdpvar(B.length,1);  V22 = BSpline(B, cV22);
cV33 = sdpvar(B.length,1);  V33 = BSpline(B, cV33);
c = [V11, V12, 0.5*bcos; V12, V22, 0.5*bsin; 0.5*bcos, 0.5*bsin, V33];
constr = set(c >= 0)
g = P2(1,1)*V11 + 2*P2(1,2)*V12 + P2(2,2)* V22 + V33;
sol = solvesdp(constr, g.integral, sdpsettings('solver', 'sdpt3'))
g1h = double(g);

figure(2), hold all
plot(p, g1h.f(p))

w = linspace(-2,2,1e2)';
figure(1), hold all
for i = 1:length(p)
    c = [bcos.f(p(i)), bsin.f(p(i))];
    nc = [-bsin.f(p(i)), bcos.f(p(i))];
    Fline = g1h.f(p(i))/norm(c)^2 * ones(size(w))*c + w*nc;
    plot(Fline(:,1), Fline(:,2), ':', 'Color', [0.6,0.6,0.6])
    plot(g1h.f(p(i))/norm(c)^2 * ones(size(w))*c(1), g1h.f(p(i))/norm(c)^2 * ones(size(w))*c(2), 'rx')
end


% Dual for radial expansion     <??? How to use splinevar here ???>
cV11 = sdpvar(B.length,1);  V11 = BSpline(B, cV11);
cV12 = sdpvar(B.length,1);  V12 = BSpline(B, cV12);
cV22 = sdpvar(B.length,1);  V22 = BSpline(B, cV22);
cV33 = sdpvar(B.length,1);  V33 = BSpline(B, cV33);
cy   = sdpvar(B.length,1);  y   = BSpline(B, cy);
c = [V11, V12, 0.5*bcos-y*bsin; V12, V22, 0.5*bsin+y*bcos; 0.5*bcos-y*bsin, 0.5*bsin+y*bcos, V33];
constr = set(c >= 0);
g = P2(1,1)*V11 + 2*P2(1,2)*V12 + P2(2,2)* V22 + V33;
sol = solvesdp(constr, g.integral, sdpsettings('solver', 'sdpt3'))
g2h = double(g);

figure(2), hold all
plot(p, g2h.f(p))

x2hd = g2h*bdir;
X2hd = cell2mat(x2hd.f(p));
X2hd = [X2hd(1:2:end), X2hd(2:2:end)];

figure(1), hold all
plot(X2hd(:,1), X2hd(:,2), '.-')