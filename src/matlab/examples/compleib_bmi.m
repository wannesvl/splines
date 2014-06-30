clear all; close all; clc;

% BMI H_inf test

% General COMPleib syntax:
name = 'REA2';
% K_range = {[1, 8], [1, 8]};  % AC7
% K_range = {[0, 10], [0, 10]};  % AC17
% K_range = {[0.5, 2], [10, 20]};  % HE1
% K_range = {[-3, -2]};  % EB2
% K_range = {[-1.5, -1]}; % NN2
K_range = {[-15, -5], [-10, -1], [0.5, 10], [-10, -0.5]};  % REA2
degree = 2;
n_knots = 2;

% Load data
[A,B1,B,C1,C,D11,D12,D21,nx,nw,nu,nz,ny] = COMPleib(name);

% Define variables
basis = @(K) BSplineBasis([K(1) * ones(1, degree), linspace(K(1), K(2), n_knots), K(2) * ones(1, degree)], degree);
bases = cellfun(basis, K_range, 'uni', false);
X = BSpline.sdpvar(bases, [nx, nx], 'symmetric');
gamma = BSpline.sdpvar(bases, [1, 1]);

% Define controller coefficients as linear functions
% Unfortunately this is not so straightforward... :-s
if nu * ny > 1
    K_coeffs = cell(2 * ones(1, nu * ny));
    [K_coeffs{:}] = deal(zeros(nu, ny));
    P = num2cell(eye(nu * ny) + 1);
    count = 1;
    for i=1:nu
        for j=1:ny
            temp = zeros(nu, ny);
            temp(i, j) = 1;
            K_coeffs{P{:, count}} = temp;
            count = count + 1;
        end
    end
else
    K_coeffs = [0, 1];
end
K = Polynomial(K_coeffs);
K = K.to_bspline(K_range);

% Define parametric optimization problem
AK = A + B * K * C;
CK = C1 + D12 * K * C;
XAK = X * AK;
XB1 = X * B1; 
LMI = [XAK' + XAK, XB1, CK';
       XB1', -gamma * eye(nw), D11';
       CK, D11, -gamma * eye(nz)];

options = sdpsettings('solver', 'sdpt3', 'sdpt3.maxit', 100, 'dualize', 0, 'verbose', 2);
sol = solvesdp([LMI <= 0, X >=0, gamma >= 0], gamma.integral, options);

% Determine local minumum in parameter range, starting point from coefficient values
gamma = double(gamma);
Jgamma = gamma.gradient;
Hgamma = gamma.hessian;
% % Maybe we could also overload max en min to get an estimate or automate this
coeffs = gamma.coeffs.coeffs2tensor;
[c_min, idx] = min(coeffs(:));
i = cell(1, nu * ny);
[i{:}] = ind2sub(size(coeffs), idx);
g = cellfun(@(b) b.greville, bases, 'uni', false);
x = cellfun(@(g, j) g(j), g, i);
options = optimoptions(@fmincon, 'GradObj','on');
obj = @(x) gamma.f(num2cell(x));
grad = @(x) Jgamma.f(num2cell(x));
hess = @(x, l) Hgamma.f(num2cell(x));
lb = cellfun(@(k) k(1), K_range);
ub = cellfun(@(k) k(2), K_range);
if nu * ny == 1
    X = fmincon(@(x) deal(obj(x), grad(x)), x, [], [], [], [], lb, ub, [], options);
else
    X = fmincon(@(x) deal(obj(x), cell2mat(grad(x))), x, [], [], [], [], lb, ub, [],  options);
end
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
    G = gamma.f({k1, k2})';
    surf(K1, K2, G)
    camlight left; light; lighting phong; shading interp;
    hold on
    plot3(X(1), X(2), gamma.f(num2cell(X)), 'r.')
    plot3(X(1), X(2), 0, 'r.')
    contour(k1, k2, G, logspace(log(min(G(:))) / log(10), log(max(G(:))) / log(10), 50))
end
