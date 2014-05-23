% Comparisons for the store3 model
%
% 1) Pólya vs knot insertion relaxation (for fixed degree of HPPD variables)
% 2) HPPD vs Bspline parameterization (increase HPPD degree vs more knots)
%
% p -> polynomial degree of the Lyapunov matrix
% n -> #internal knots spline LMI variables
% m -> #midpoint refinements, s.t. #inserted knots = 2^m - 1
% d -> degree of Pólya relaxation

clear all
close all
clear global
clc

% [Prim,Dual] = store3_splinevars(1,0,0,0)

% 1) Pólya vs knot insertion 
p = 3; n = 0; m = 0; % Polya
for d = 0:4 
    [Prim,Dual]   = store3_splinevars(p,n,m,2^d-1);
    Polya.p(d+1)  = Prim.obj;
    Polya.cp(d+1) = Prim.ctime;
    Polya.d(d+1)  = Dual.obj;
    Polya.cd(d+1) = Dual.ctime;
end
p = 3; n = 0; d = 0; % knot insertion
count = 1;
for m = 0:4 
    [Prim,Dual]  = store3_splinevars(p,n,m,d);
    KI.p(count)  = Prim.obj;
    KI.cp(count) = Prim.ctime;
    KI.d(count)  = Dual.obj;
    KI.cd(count) = Dual.ctime;
    count = count + 1;
end

% 2) HPPD vs Bspline parameterization
m = 0; d = 0; n = 0; % increase polynomial degree
for p = 2:5 
    [Prim,Dual] = store3_splinevars(p,n,m,d);
    deg.p(p-1)  = Prim.obj;
    deg.cp(p-1) = Prim.ctime;
    deg.d(p-1)  = Dual.obj;
    deg.cd(p-1) = Dual.ctime;
end
m = 0; d = 0; p = 2; % increase #knots
for n = 0:3
    [Prim,Dual]   = store3_splinevars(p,n,m,d);
    knots.p(n+1)  = Prim.obj;
    knots.cp(n+1) = Prim.ctime;
    knots.d(n+1)  = Dual.obj;
    knots.cd(n+1) = Dual.ctime;
end
