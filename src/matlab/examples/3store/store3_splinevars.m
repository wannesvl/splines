clear all
% close all
clear global
clc

%_______________________________________________________________________
% STEP 1. Load model and sampled solution
load model_3store2
load gamma_k2
mu = 15;

% conversion system A matrix to Bspline
temp1 = A{2,1}; temp2 = A{1,2};
clear A
cA = {temp2;temp1}; % Bspline coefficients of A
A  = BSpline(BSplineBasis([0 0 1 1],1),cA); 

%_______________________________________________________________________      
% STEP 2: solve parametric LMI optimization problems
p = 4; % degree of Lyapunov matrix Q
n = 5; % #internal knots of LMI variables Q,Z,L
m = 0; % #inserted knots 
d = 0; % degree Pólya relaxation

tol = 1e-5; % tolerance primal residual
eps = 0;    % tuning parameter to enforce strict inequality LMIs

%_______________________________________________________________________ 
% 2.1: Solve parametric primal

% new LMI system
LMIs = set([]);
Primal.rows = 0;

% generate the LMI variables
BasisQ  = BSplineBasis([0*ones(1,p)  , linspace(0,1,n+2), ones(1,p)  ],p);
BasisZL = BSplineBasis([0*ones(1,p+1), linspace(0,1,n+2), ones(1,p+1)],p+1);
for j = 1:length(BasisQ)
    cQ{j,:} = sdpvar(nx,nx,'symmetric');
end
Q = BSpline(BasisQ,cQ);
for j = 1:length(BasisZL)
    cZ{j,:} = sdpvar(nu,nu,'symmetric');
    cL{j,:} = sdpvar(nu,nx,'full');
end
Z = BSpline(BasisZL,cZ);
L = BSpline(BasisZL,cL);

% construct the LMI terms
Term1 = [Q*A' + A*Q + Bu*L + L'*Bu', Bw; Bw', -eye(nw)]; 
Term2 = Cz*Q*Cz' - mu*eye(nz);
Term3 = [Q, L'; L, Z];

% Pólya relaxations
Term1 = Term1.increase_degree(d);
Term2 = Term2.increase_degree(1+d);
Term3 = Term3.increase_degree(d);

% knot insertion
insknots = linspace(0,1,m+2); insknots(1) = []; insknots(end) = [];
Term1 = Term1.insert_knots(insknots);
Term2 = Term2.insert_knots(insknots);
Term3 = Term3.insert_knots(insknots);

% sufficient LMIs (all Bspline coefficients >/< 0)
for j = 1:length(Term1.basis)
    LMIs = LMIs + set(Term1.getcoeffs{j} < -eps*eye(nx+nw));
    Primal.rows = Primal.rows + nx + nw;
end
for j = 1:length(Term2.basis) 
    LMIs = LMIs + set(Term2.getcoeffs{j} < -eps*eye(nz));
    Primal.rows = Primal.rows + nz;
end
for j = 1:length(Term3.basis) 
    LMIs = LMIs + set(Term3.getcoeffs{j} > eps*eye(nx + nu));
    Primal.rows = Primal.rows + nx + nu;
end
Primal.vars = size(getvariables(LMIs),2); % #LMI variables

% objective function
obj = Z.trace.integral;

% solve the SDP
sol = solvesdp(LMIs, obj, sdpsettings('solver', 'sdpt3', 'sdpt3.maxit', 100, 'dualize', 0, 'verbose', 1))
[p_res,d_res] = checkset(LMIs);

% extract solution variables
Primal.feas = 0;
if p_res > -tol
    Primal.feas = 1;
    cQ = cellfun(@double, cQ, 'uni', false);
    cZ = cellfun(@double, cZ, 'uni', false);
    cL = cellfun(@double, cL, 'uni', false);
    Primal.obj = double(obj);
end

%_______________________________________________________________________
% 2.2: Solve parametric dual

% new LMI system
LMIs = set([]);
Dual.rows = 0;

% generate the LMI variables
BasisUV = BSplineBasis([0*ones(1,p), linspace(0,1,n+2), ones(1,p)],p);
for j = 1:length(BasisUV)
    cU11{j,:} = sdpvar(nx,nx,'symmetric');
    cU12{j,:} = sdpvar(nx,nw,'full');
    cU22{j,:} = sdpvar(nw,nw,'symmetric');
    cV{j,:} = sdpvar(nz,nz,'symmetric');
end
U11 = BSpline(BasisUV,cU11);
U12 = BSpline(BasisUV,cU12);
U22 = BSpline(BasisUV,cU22);
U = [U11, U12; U12', U22];
V = BSpline(BasisUV,cV);

% construct the LMI terms
Term = [A'*U11 + U11*A + Cz'*V*Cz, U11*Bu; Bu'*U11, eye(nu)]; 

% Pólya relaxation
Term = Term.increase_degree(d);

% knot insertion
insknots = linspace(0,1,m+2); insknots(1) = []; insknots(end) = [];
Term = Term.insert_knots(insknots);

% sufficient LMIs (all Bspline coefficients >/< 0)
for j = 1:length(U11.basis)
    LMIs = LMIs + set(U.getcoeffs{j} > eps*eye(nx+nw));
    Dual.rows = Dual.rows + nx + nw;
    LMIs = LMIs + set(V.getcoeffs{j} > eps*eye(nz));
    Dual.rows = Dual.rows + nz;
end
for j = 1:length(Term.basis)
    LMIs = LMIs + set(Term.getcoeffs{j} > eps*eye(nx+nu));
    Dual.rows = Dual.rows + nx + nu;
end
Dual.vars = size(getvariables(LMIs),2); % #LMI variables

% objective function
obj1 = U12*Bw'; obj2 = U12'*Bw; obj3 = -U22; obj4 = -mu*V;
obj = obj1.trace.integral + obj2.trace.integral + obj3.trace.integral + obj4.trace.integral;

% solve the SDP
sol = solvesdp(LMIs, -obj, sdpsettings('solver', 'sdpt3', 'sdpt3.maxit', 100, 'dualize', 0, 'verbose', 1))
[p_res,d_res] = checkset(LMIs);

% extract solution variables
Dual.feas = 0;
if p_res > -tol
    Dual.feas = 1;
    cU11 = cellfun(@double, cU11, 'uni', false);
    cU12 = cellfun(@double, cU12, 'uni', false);
    cU22 = cellfun(@double, cU22, 'uni', false);
    cV   = cellfun(@double, cV, 'uni', false);
    Dual.obj = double(obj);
end

% objective function vs parameter
U11 = BSpline(BasisUV,cU11);
U12 = BSpline(BasisUV,cU12);
U22 = BSpline(BasisUV,cU22);
V = BSpline(BasisUV,cV);
obj1 = U12*Bw'; obj2 = U12'*Bw; obj3 = -U22; obj4 = -mu*V;
obj = obj1.trace + obj2.trace + obj3.trace + obj4.trace;

%_______________________________________________________________________
% 2.3: Plot results

% sampled solution
figure(1); plot(k2_vec, gamma_vec); hold on;

% primal objective vs parameter
Z = BSpline(BasisZL, cZ); % redefine B-spline in terms of double coeffs
theta = linspace(0, 1, 101);
figure(1); plot(k2L+(k2U-k2L)*theta, Z.trace.f(theta),'r','Linewidth',1.5)

% dual objective vs parameter
figure(1); plot(k2L+(k2U-k2L)*theta, obj.f(theta),'g','Linewidth',1.5)

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

