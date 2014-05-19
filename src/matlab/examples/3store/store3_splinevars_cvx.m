clear all
% close all
clc

%_______________________________________________________________________
% STEP 1. Load model and sampled solution
load model_3store2
load gamma_k2
mu = 15;

% conversion system matrices to Bsplines
temp1 = A{2,1}; temp2 = A{1,2};
clear A
cA = {temp2; temp1}; % Bspline coefficients of A
A  = BSpline(BSplineBasis([0 0 1 1], 1), cA);

%_______________________________________________________________________      
% STEP 2: solve parametric LMI optimization problems

p = 2; % degree of Lyapunov matrix Q
% n = 0; % #knots of Q
BasisQ  = BSplineBasis([0 * ones(1, p + 1) ones(1, p + 1)], p);
BasisZL = BSplineBasis([0 * ones(1, p + 2) ones(1, p + 2)], p + 1);

%_______________________________________________________________________ 
% 2.1: Solve parametric primal
cvx_solver sedumi
cvx_begin sdp
    variable cZ(nu, nu, length(BasisZL)) symmetric;
    variable cL(nu, nx, length(BasisZL));
    variable cQ(nx, nx, length(BasisQ)) symmetric;
        
    Q = BSpline(BasisQ, squeeze(mat2cell(cQ, nx, nx, ones(length(BasisQ), 1))));
    Z = BSpline(BasisZL, squeeze(mat2cell(cZ, nu, nu, ones(length(BasisZL), 1))));
    L = BSpline(BasisZL, squeeze(mat2cell(cL, nu, nx, ones(length(BasisZL), 1))));

    LMI1 = [Q * A' + A * Q + Bu * L + L' * Bu', Bw; Bw', -1]; % LMI1 < 0
    LMI2 = Cz * Q * Cz' - mu * eye(nz); % LMI2 < 0
    LMI2 = LMI2.increase_degree(1)
    LMI3 = [Q, L'; L, Z]; % LMI3 > 0

    minimize (Z.trace.integral)
    subject to
        for i=1:length(LMI1.basis)
            LMI1.coeffs.coeffs{i} < 0
        end
        for i=1:length(LMI2.basis)
            LMI2.coeffs.coeffs{i} < 0
        end
        for i=1:length(LMI3.basis)
            LMI3.coeffs.coeffs{i} > 0
        end
cvx_end

Z = BSpline(BasisZL, squeeze(mat2cell(cZ, nu, nu, ones(length(BasisZL), 1))));
theta = linspace(0, 1, 101);

figure(1)
hold on
plot(theta, Z.trace.f(theta))
