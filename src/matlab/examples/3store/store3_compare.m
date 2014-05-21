% Comparisons for the store3 model
% 1) Pólya vs knot insertion relaxation (for fixed degree of HPPD variables)
% 2) HPPD vs Bspline parameterization (increase HPPD degree vs more knots)

clear all
close all
clear global
clc

% [Prim,Dual] = store3_splinevars(1,0,1,1);

% 1) Pólya vs knot insertion 
p = 2; % degree of Lyapunov matrix Q
n = 0; % #internal knots of LMI variables Q,Z,L
for m = 0:3 % #inserted knots
    for d = 0:3 % degree of Pólya relaxation
        [Prim,Dual] = store3_splinevars(p,n,m,d);
        g1_prim(m+1,d+1) = Prim.obj;
        g1_dual(m+1,d+1) = Dual.obj;
    end
end

% 2) HPPD vs Bspline parameterization
m = 0; % #inserted knots
d = 0; % degree of Pólya relaxation 
for p = 1:4 % degree of Lyapunov matrix Q
    for n = 0:3 % #internal knots of LMI variables Q,Z,L
        [Prim,Dual] = store3_splinevars(p,n,m,d);
        g2_prim(p,n+1) = Prim.obj;
        g2_dual(p,n+1) = Dual.obj;
    end
end