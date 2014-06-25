clear all; close all; clc;

% BMI H_inf test

% General COMPleib syntax:
name = 'AC7';
K_range = {[1, 8], [1, 8]};  % AC7
% K_range = {[-3, -2]};  % EB2
% K_range = {[-1.5, -1]}; % NN2
degree = 2;
n_knots = 5;

% Load data
[A,B1,B,C1,C,D11,D12,D21,nx,nw,nu,nz,ny] = COMPleib(name);

% Define variables
basis = @(K) BSplineBasis([K(1) * ones(1, degree), linspace(K(1), K(2), n_knots), K(2) * ones(1, degree)], degree);
bases = cellfun(basis, K_range, 'uni', false);
X = BSpline.sdpvar(bases, [nx, nx], 'symmetric');
gamma = BSpline.sdpvar(bases, [1, 1]);

% Define controller coefficients as linear functions
% Sadly this is not so straightforward... :-s
if nu * ny > 1
    K_coeffs = cell(2 * ones(1, nu * ny));
    [K_coeffs{:}] = deal(zeros(nu, ny));
    P = eye(nu * ny);
    P = num2cell(P + 1);
    for i=1:nu * ny
        temp = zeros(1, nu * ny);
        temp(i) = 1;
        K_coeffs{P{:, i}} = temp;
    end
else
    K_coeffs = [0, 1];
end
K = Polynomial(K_coeffs);
K = K.to_bspline(K_range);

% Define parametric optimization problem
AK = A + B * K * C;
CK = C1 + D12 * K * C;
LMI = [AK' * X + X * AK, X * B1, CK';
       B1' * X, -gamma * eye(nw), D11';
       CK, D11, -gamma * eye(nz)];

options = sdpsettings('solver', 'sdpt3', 'sdpt3.maxit', 100, 'dualize', 0, 'verbose', 2);
sol = solvesdp([LMI <= 0, X >=0, gamma >= 0], gamma.integral, options);

% Determine local minumum in parameter range, starting point from coefficient values
gamma = double(gamma);
% Maybe we could also overload max en min to get an estimate or automate this
coeffs = gamma.coeffs.coeffs2tensor;
[c_min, idx] = min(coeffs(:));
i = cell(1, nu * ny);
[i{:}] = ind2sub(size(coeffs), idx);
g = cellfun(@(b) b.greville, bases, 'uni', false);
x = cellfun(@(g, j) g(j), g, i);
X = fminunc(@(x) gamma.f(num2cell(x)), x);
X, gamma.f(num2cell(X))

% Plot for 1D or 2D problems
if nu * ny == 1
    figure
    k = K_range{1};
    k = linspace(k(1), k(2), 101);
    plot(k, gamma.f(k)); hold on;
    plot(X, gamma.f(X), 'rx')
elseif nu * ny == 2
    figure
    k1 = K_range{1};
    k2 = K_range{2};
    k1 = linspace(k1(1), k1(2), 101);
    k2 = linspace(k2(1), k2(2), 101);
    [K1, K2] = meshgrid(k1, k2);
    surf(K1, K2, gamma.f({k1, k2}))
    camlight left; light; lighting phong; shading interp;
    hold on
    plot3(X(1), X(2), gamma.f(num2cell(X)))
end
