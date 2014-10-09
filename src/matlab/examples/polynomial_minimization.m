close all; clear all; clc;

% Motzkin polynomial
f = @(x,y) x.^2 .* y.^2 .* (x.^2 + y.^2 - 1) + 1/27;
figure
[X,Y] = meshgrid([-1:0.02:1],[-1:0.02:1]);
contour(X,Y,f(X,Y),100)

% Yalmip moment relaxations
x = sdpvar(1, 1);
y = sdpvar(1, 1);
h = x^2 * y^2 * (x^2 + y^2 - 1) + 1/27;
F = [];

moment = [];
moment_time = [];
for m=0:5
    for i=1:5
        sol = solvemoment(F, h, sdpsettings('solver','sdpt3'), m+3);
        time(i) = sol.solvertime;
    end
    moment(m+1) = relaxdouble(h);
    moment_time(m+1) = mean(time);
end

% BSpline
% =======
d = 3;
A = zeros(4, 4);
A(1,1) = 1/27; A(3,3) = -1; A(5,3) = 1; A(3,5) = 1;
P = Polynomial(A);
P = P.to_bspline({[-1.5,1.5],[-1.5,1.5]});

% Degree elevation polya
t = sdpvar(1, 1);
polya = [];
polya_time = [];
for i=0:9
    Pi = P.increase_degree([i, i]);
    con = [Pi-t>=0];
    for j=1:10
        sol = solvesdp(con, -t, sdpsettings('solver','clp'));
        time(j) = sol.solvertime;
    end
    polya(i+1) = double(t);
    polya_time(i+1) = mean(time);
end

% Knot insertion
bspline = [polya(1)];
bspline_time = [polya_time(1)];
for i=1:6
    Pi = P.insert_knots({linspace(-1.5, 1.5, 2^i+1), linspace(-1.5, 1.5, 2^i+1)});
    con = [Pi-t>=0];
    for j=1:10
        sol = solvesdp(con, -t, sdpsettings('solver','clp'));
        time(j) = sol.solvertime;
    end
    bspline(i+1) = double(t);
    bspline_time(i+1) = mean(time);
end

figure
semilogx(moment_time, moment)
hold all
semilogx(polya_time, polya)
semilogx(bspline_time, bspline)

figure
loglog(moment_time, -moment)
hold all
loglog(polya_time, -polya)
loglog(bspline_time, -bspline)