function [Primal,Dual] = store3_splines_mvar(p,n,m,d)
%
% Combined structure and control design for the 3 store problem.
%
% State feedback synthesis.
% Structural parameters: k2 (spring) and d3 (damping)
%
% The objective is to solve only one primal and one dual optimization 
% problem to obtain a good approximation of the closed-loop performance as 
% a function of the parameter (k2,d3). 
%
% This problem is formulated as an LPV state-feedback synthesis problem 
% with constant and measurable parameters, giving rise to a semi-infinite 
% convex optimization problem with LMI constraints.
%
% The LMI variables are assumed to have a multivariate B-spline parameter 
% dependency.
% 
% Inputs:   It is assumed that the following are identical for each basis:
%           p -> polynomial degree of the Lyapunov matrix 
%           n -> #internal knots spline LMI variables 
%           m -> #midpoint refinements, s.t. #inserted knots = 2^m - 1
%           d -> degree of Polya relaxation
%
% Outputs:  Primal: struct with info about primal solution
%           Dual  : struct with info about dual solution
%
%
% References: 
%
% J. F. Camino, M. C. de Oliveira and R. E. Skelton. "Convexifying" Linear 
% Matrix Inequality Methods for Integrating Structure and Control Design. 
% Journal of Structural Engineering, July 2003.
%
%_______________________________________________________________________
% STEP 1: set problem parameters 
mu = 15;

%_______________________________________________________________________
% STEP 2: Sampled solution as a function of (k2,d3)
% load model_3store2_mvar
% 
% k2_vec = linspace(k2L, k2U, 30)';
% d3_vec = linspace(d3L, d3U, 30)';
% gamma_mat = zeros(length(k2_vec),length(d3_vec));
% 
% count = 0;
% for i = 1:length(k2_vec)
%     for j = 1:length(d3_vec)
%         k2 = k2_vec(i);
%         d3 = d3_vec(j);
%         
%         K = [k1+k2, -k2  , 0  ;
%              -k2  , k2+k3, -k3;
%              0    , -k3  , k3 ];
%         D = [d1+d2, -d2  , 0  ;
%              -d2  , d2+d3, -d3;
%              0    , -d3  , d3 ];      
%         A = Ei*[zeros(3), eye(3); -K, -D];
% 
%         Q = sdpvar(nx); L = sdpvar(nu,nx); Z = sdpvar(nu);
%         eps = 1e-9;
% 
%         constr = [Q >= eps*eye(nx), [Q*A'+A*Q'+Bu*L+L'*Bu', Bw; Bw', -eye(nw)] <= 0];
%         constr = [constr, Cz*Q*Cz' <= mu*eye(nz)];
%         constr = [constr, [Q, L'; L, Z] >= 0];
%         sol = solvesdp(constr, trace(Z), sdpsettings('solver', 'sdpt3', 'sdpt3.maxit', 100, 'dualize', 1, 'verbose', 0));
%         gamma_mat(i,j) = double(trace(Z));
%         count = count + 1
%     end
% end

% save('sampled_solution_multivar.mat','k2_vec','d3_vec','gamma_mat');

%_______________________________________________________________________      
% STEP 3: solve parametric LMI optimization problems
load model_3store2_mvar % reload B-spline model
tol = 1e-5;             % tolerance primal residual
eps = 0;                % tuning parameter strict inequality LMIs

%_______________________________________________________________________ 
% 3.1: Solve parametric primal
tstart = tic;

% define bases for LMI variables
BasisQk = BSplineBasis([k2L*ones(1,p), linspace(k2L,k2U,n+2), k2U*ones(1,p)], p);
BasisQd = BSplineBasis([d3L*ones(1,p), linspace(d3L,d3U,n+2), d3U*ones(1,p)], p);
BasisZk = BSplineBasis([k2L*ones(1,p), linspace(k2L,k2U,n+2), k2U*ones(1,p)], p);
BasisZd = BSplineBasis([d3L*ones(1,p), linspace(d3L,d3U,n+2), d3U*ones(1,p)], p);
BasisLk = BSplineBasis([k2L*ones(1,p), linspace(k2L,k2U,n+2), k2U*ones(1,p)], p);
BasisLd = BSplineBasis([d3L*ones(1,p), linspace(d3L,d3U,n+2), d3U*ones(1,p)], p);

% generate the LMI variables
Q = BSpline.sdpvar({BasisQk,BasisQd}, [nx, nx], 'symmetric');
Z = BSpline.sdpvar({BasisZk,BasisZd}, [nu, nu], 'symmetric');
L = BSpline.sdpvar({BasisLk,BasisLd}, [nu, nx], 'full');

% construct the LMI terms
Term1 = [Q*A' + A*Q' + Bu*L + L'*Bu', Bw; Bw', -eye(nw)]; 
Term2 = Cz*Q*Cz' - mu*eye(nz);
Term3 = [Q, L'; L, Z];

% knot insertion
insknots_k2 = linspace(k2L,k2U,2^m+1); insknots_k2(1) = []; insknots_k2(end) = [];
insknots_d3 = linspace(d3L,d3U,2^m+1); insknots_d3(1) = []; insknots_d3(end) = [];
Term1 = Term1.insert_knots({insknots_k2,insknots_d3});
Term2 = Term2.insert_knots({insknots_k2,insknots_d3});
Term3 = Term3.insert_knots({insknots_k2,insknots_d3});

% Polya relaxations: increase degree of each basis by d
Term1 = Term1.increase_degree([d,d]);
Term2 = Term2.increase_degree([d,d]);
Term3 = Term3.increase_degree([d,d]);

% solve the SDP
obj     = Z.trace.integral;
con     = [Term1 <= 0, Term2 <= 0, Term3 >= 0];
options = sdpsettings('solver', 'sdpt3', 'sdpt3.maxit', 100, 'dualize', 0, 'verbose', 2);
sol     = solvesdp(con, obj, options);

% Get SDP info
info_        = lmiinfo(con);
Primal.rows  = sum(info_.sdp(:, 1));
Primal.vars  = size(getvariables(con), 2);
Primal.ctime = toc(tstart);
[pr,~] = checkset(con); % primal residual
if pr > -tol
    Primal.obj = double(obj);
else
    Primal.obj = [];
    return;
end

%_______________________________________________________________________ 
% 3.2: Solve parametric dual
tstart = tic;

% define bases for LMI variables
BasisUVk = BSplineBasis([k2L*ones(1,p), linspace(k2L,k2U,n+2), k2U*ones(1,p)], p);
BasisUVd = BSplineBasis([d3L*ones(1,p), linspace(d3L,d3U,n+2), d3U*ones(1,p)], p);

% generate the LMI variables
U11 = BSpline.sdpvar({BasisUVk,BasisUVd}, [nx, nx], 'symmetric');
U12 = BSpline.sdpvar({BasisUVk,BasisUVd}, [nx, nw], 'full');
U22 = BSpline.sdpvar({BasisUVk,BasisUVd}, [nw, nw], 'symmetric');
U = [U11, U12; U12', U22];
V = BSpline.sdpvar({BasisUVk,BasisUVd}, [nz, nz], 'symmetric');

% construct the LMI terms
Term = [A'*U11 + U11*A + Cz'*V*Cz, U11*Bu; Bu'*U11, eye(nu)]; 

% knot insertion
insknots_k2 = linspace(k2L,k2U,2^m+1); insknots_k2(1) = []; insknots_k2(end) = [];
insknots_d3 = linspace(d3L,d3U,2^m+1); insknots_d3(1) = []; insknots_d3(end) = [];
Term = Term.insert_knots({insknots_k2,insknots_d3});

% Polya relaxation: increase degree of each basis by d
Term = Term.increase_degree([d,d]);

% solve the SDP
obj1    = U12*Bw'; obj2 = U12'*Bw; obj3 = -U22; obj4 = -mu*V;
obj     = obj1.trace + obj2.trace + obj3.trace + obj4.trace;
obj     = obj.integral;
con     = [U >= 0, V >= 0, Term >= 0];
options = sdpsettings('solver', 'sdpt3', 'sdpt3.maxit', 100, 'dualize', 0, 'verbose', 2);
sol     = solvesdp(con, -1e-15*obj, options);

% Get SDP info
info_      = lmiinfo(con);
Dual.rows  = sum(info_.sdp(:, 1));
Dual.vars  = size(getvariables(con), 2);
Dual.ctime = toc(tstart);
[pr,~] = checkset(con); % primal residual
if pr > -tol
    Dual.obj = double(obj);
else
    Dual.obj = [];
    return;
end

%_______________________________________________________________________      
% STEP 4: generate plots

% sampled solution
load sampled_solution_mvar.mat
[x,y] = meshgrid(k2_vec, d3_vec);
[gamma_min,ind_min] = min(gamma_mat(:));          % optimal objective
[ind_x,ind_y] = ind2sub(size(gamma_mat),ind_min); % corresp. index (x,y)
x_min = x(1,ind_x); y_min = y(ind_y,1);           % corresp. (x,y)
figure(1); surf(x,y,gamma_mat'); hold on; plot3(x_min,y_min,gamma_min,'g*');

% primal objective vs parameter
dxy = 1/50; Z = double(Z);
[x,y] = meshgrid(k2L:(k2U-k2L)*dxy:k2U, d3L:(d3U-d3L)*dxy:d3U);
z = Z.trace.f({x(1,:),y(:,1)})';
[z_min,ind_min] = min(z(:));              % minimum objective
[ind_y,ind_x] = ind2sub(size(z),ind_min); % corresponding index (x,y)
x_min = x(1,ind_x); y_min = y(ind_y,1);   % corresponding (x,y)
figure(1); surf(x,y,z); hold on; plot3(x_min,y_min,z_min,'r*');
xlabel('k2'); ylabel('d3'); zlabel('H2 Bound'); title('Primal & Dual solution');

% dual objective vs parameter
obj1 = double(obj1); obj2 = double(obj2); obj3 = double(obj3); obj4 = double(obj4);
obj = obj1.trace + obj2.trace + obj3.trace + obj4.trace;
z = obj.f({x(1,:),y(:,1)})';
figure(1); surf(x,y,z);

