clear all; close all; clc;

% BMI H_inf test

% General COMPleib syntax:
name = 'HE1';
% K_range = {[1, 9], [1, 8]};  % AC7
% K_range = {[0, 10], [0, 10]};  % AC17
K_range = {[0.5, 2], [10, 20]};  % HE1
% K_range = {[-3, -2]};  % EB2
% K_range = {[-1.5, -1]}; % NN2
% K_range = {[-15, -5], [-10, -1], [0.5, 10], [-10, -0.5]};  % REA2
degree = 2;
n_knots = 8;

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
[K, m] = min(gamma);
K = reshape(K, nu, ny);

% Evaluate true gamma for chosen control parameters
X = sdpvar(nx, nx, 'symmetric');
g = sdpvar(1, 1);
AK = A + B * K * C;
CK = C1 + D12 * K * C;
XAK = X * AK;
XB1 = X * B1; 
LMI = [XAK' + XAK, XB1, CK';
       XB1', -g * eye(nw), D11';
       CK, D11, -g * eye(nz)];

options = sdpsettings('solver', 'sdpt3', 'sdpt3.maxit', 100, 'dualize', 0, 'verbose', 2);
sol = solvesdp([LMI <= 0, X >=0, g >= 0], g, options);

% Show results to the user
mes = sprintf('gamma: %0.5g \nreevaluation: %0.5g', m, double(g));
disp(mes)

% Plot for 1D or 2D problems
if nu * ny == 1
    figure
    k = K_range{1};
    k = linspace(k(1), k(2), 101);
    plot(k, gamma.f(k)); hold on;
    plot(K, m, 'rx')
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
    plot3(K(1), K(2), m, 'r.')
    plot3(K(1), K(2), 0, 'r.')
    contour(k1, k2, G, logspace(log(min(G(:))) / log(10), log(max(G(:))) / log(10), 50))
end
