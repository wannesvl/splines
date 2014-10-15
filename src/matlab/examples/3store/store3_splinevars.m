function [Primal,Dual] = store3_splinevars(p,n,m,d)
% State-feedback synthesis for the 3 store problem.
% The spring coefficient k2 is a structural parameter, 
% and is assumed to be constant and measurable.
%
% The objective is to solve only one primal and one dual optimization 
% problem to obtain a good approximation of the performance as a function
% of the parameter. Therefore, this problem is formulated as an LPV
% state-feedback synthesis problem with a constant and measurable
% parameter, giving rise to a semi-infinite convex optimization problem 
% with LMI constraints. 
%
% The LMI variables are assumed 
% 
% Inputs:   p -> polynomial degree of the Lyapunov matrix
%           n -> #internal knots spline LMI variables
%           m -> #midpoint refinements, s.t. #inserted knots = 2^m - 1
%           d -> degree of Pólya relaxation
%
% Outputs:  Primal: struct with info about primal solution
%           Dual  : struct with info about dual solution
%
%_______________________________________________________________________
% STEP 1. Load model and sampled solution

load model_3store2
load gamma_k2
mu = 15;

% conversion system A matrix to Bspline
A  = BSpline(BSplineBasis([0, 0, 1, 1], 1), {A{1, 2}, A{2, 1}}); 
% A = Polynomial({A{1, 2}, A{2, 1} - A{1, 2}});  % This does not seem to work...

%_______________________________________________________________________      
% STEP 2: solve parametric LMI optimization problems
tol = 1e-5; % tolerance primal residual
eps = 0;    % tuning parameter to enforce strict inequality LMIs

%_______________________________________________________________________ 
% 2.1: Solve parametric primal
tstart = tic;

BasisQ  = BSplineBasis([0*ones(1,p)  , linspace(0,1,n+2), ones(1,p)  ],p);
BasisZL = BSplineBasis([0*ones(1,p+1), linspace(0,1,n+2), ones(1,p+1)],p+1);
% generate the LMI variables
Q = BSpline.sdpvar(BasisQ, [nx, nx], 'symmetric');
Z = BSpline.sdpvar(BasisZL, [nu, nu], 'symmetric');
L = BSpline.sdpvar(BasisZL, [nu, nx], 'full');

% construct the LMI terms
Term1 = [Q*A' + A*Q' + Bu*L + L'*Bu', Bw; Bw', -eye(nw)]; 
Term2 = Cz*Q*Cz' - mu*eye(nz);
Term3 = [Q, L'; L, Z];

% knot insertion
insknots = linspace(0,1,2^m+1); insknots(1) = []; insknots(end) = [];
Term1 = Term1.insert_knots({insknots});
Term2 = Term2.insert_knots({insknots});
Term3 = Term3.insert_knots({insknots});

% % Pólya relaxations
Term1 = Term1.increase_degree(d);
Term2 = Term2.increase_degree(d);
Term3 = Term3.increase_degree(d);

% objective function
obj = Z.trace.integral;

% solve the SDP
con = [Term1 <= 0, Term2 <= 0, Term3 >= 0];
options = sdpsettings('solver', 'sdpt3', 'sdpt3.maxit', 100, 'dualize', 0, 'verbose', 2);
sol = solvesdp(con, obj, options);

% Get SDP info
info_ = lmiinfo(con);
Primal.rows = sum(info_.sdp(:, 1));
Primal.vars = size(getvariables(con), 2);
Primal.obj = double(obj);
Primal.ctime = toc(tstart);

%_______________________________________________________________________ 
% 2.2: Solve parametric primal

tstart = tic;

% generate the LMI variables
BasisUV = BSplineBasis([0*ones(1,p), linspace(0,1,n+2), ones(1,p)],p);
U11 = BSpline.sdpvar(BasisUV, [nx, nx], 'symmetric');
U12 = BSpline.sdpvar(BasisUV, [nx, nw], 'full');
U22 = BSpline.sdpvar(BasisUV, [nw, nw], 'symmetric');
U = [U11, U12; U12', U22];
V = BSpline.sdpvar(BasisUV, [nz, nz], 'symmetric');

% construct the LMI terms
Term = [A'*U11 + U11*A + Cz'*V*Cz, U11*Bu; Bu'*U11, eye(nu)]; 

% Pólya relaxation
Term = Term.increase_degree(d);

% knot insertion
insknots = linspace(0,1,2^m+1); insknots(1) = []; insknots(end) = [];
Term = Term.insert_knots({insknots});

% sufficient LMIs (all Bspline coefficients >/< 0)
con = [U >= 0, V >=0, Term >= 0];

% objective function
obj1 = U12*Bw'; obj2 = U12'*Bw; obj3 = -U22; obj4 = -mu*V;
obj = obj1.trace.integral + obj2.trace.integral + obj3.trace.integral + obj4.trace.integral;

% solve the SDP
sol = solvesdp(con, -obj, options);

% Get SDP info
info_ = lmiinfo(con);
Dual.rows = sum(info_.sdp(:, 1));
Dual.vars = size(getvariables(con), 2);
Dual.obj = double(obj);
Dual.ctime = toc(tstart);

%_______________________________________________________________________
% 2.3: Plot results

% sampled solution
figure(1); plot(k2_vec, gamma_vec); hold on;

% primal objective vs parameter
Z = double(Z); % redefine B-spline in terms of double coeffs
theta = linspace(0, 1, 101);
figure(1); plot(k2L + (k2U - k2L) * theta, Z.trace.f(theta),'r','Linewidth',1.5)

% dual objective vs parameter
obj1 = double(obj1); obj2 = double(obj2); obj3 = double(obj3); obj4 = double(obj4);
obj = obj1.trace + obj2.trace + obj3.trace + obj4.trace;
figure(1); plot(k2L + (k2U - k2L) * theta, obj.f(theta),'g','Linewidth',1.5)

% upper/lower bounds of parameter and knots
yL = get(gca,'YLim');
line([k2L k2L],yL,'color','black');
line([k2U k2U],yL,'color','black');
intknots = linspace(0,1,n+2); intknots(1) = []; intknots(end) = [];
for j = 1:length(intknots) % internal knots LMI vars
   line([k2L+(k2U-k2L)*intknots(j), k2L+(k2U-k2L)*intknots(j)],yL,'color','red');  
end
for j = 1:length(insknots) % inserted knots
   line([k2L+(k2U-k2L)*insknots(j), k2L+(k2U-k2L)*insknots(j)],yL,'color','black');  
end
