function [Primal,Dual] = store3_HPPDclass(p,d)

%_______________________________________________________________________
% STEP 1. Load model and sampled solution
load model_3store2.mat;
load gamma_k2
mu = 15;

%_______________________________________________________________________      
% STEP 2: solve parametric LMI optimization problems
tol = 1e-5;             % tolerance primal residual
A   = HPPDmatrix(A);    % system matrices
Bw  = HPPDmatrix(Bw); 
Bu  = HPPDmatrix(Bu); 
Cz  = HPPDmatrix(Cz);

%_______________________________________________________________________ 
% 2.1: Solve parametric primal

% new LMI system
LMIs = set([]);

% generate the LMI variables
KNp = get_KNg(A.vertices,p); 
for j=1:get_JNg(A.vertices,p)
    jj = get_monomial_string(KNp(j,:));
    eval(['cQ{' jj '} = sdpvar(nx,nx,''symmetric'');']);
end
Q = HPPDmatrix(cQ);
KNp1 = get_KNg(A.vertices,p+1);  
for j=1:get_JNg(A.vertices,p+1)
    jj = get_monomial_string(KNp1(j,:));
    eval(['cZ{' jj '} = sdpvar(nu,nu,''symmetric'');']);
    eval(['cL{' jj '} = sdpvar(nu,nx,''full'');']);
end
Z  = HPPDmatrix(cZ);
L  = HPPDmatrix(cL);
Iw = HPPDmatrix(eye(nw));
Iz = HPPDmatrix(eye(nz));
mu = HPPDmatrix(mu);

% construct the LMI terms
U = [Q*A' + A*Q + Bu*L + L'*Bu', Bw; Bw', -Iw];
V = Cz*Q*Cz' - mu*Iz;
V = V.increase_degree(1);
W = [Q, L'; L, Z];

% apply Pólya relaxation of degree d
U = U.increase_degree(d);
V = V.increase_degree(d);
W = W.increase_degree(d);

% sufficient LMIs
for j = 1:length(Q.getcoeffs);
    LMIs = LMIs + set(Q.getcoeffs{j} > 0);
end
for j = 1:length(U.getcoeffs)
    LMIs = LMIs + set(U.getcoeffs{j} < 0);
end
for j = 1:length(V.getcoeffs)
    LMIs = LMIs + set(V.getcoeffs{j} < 0);
end
for j = 1:length(W.getcoeffs)
    LMIs = LMIs + set(W.getcoeffs{j} > 0);
end

% objective function
goal = Z.trace.integral;

% solve the SDP
sol = solvesdp(LMIs, goal, sdpsettings('solver', 'sdpt3', 'sdpt3.maxit', 100, 'dualize', 0, 'verbose', 1));
[p_res,d_res] = checkset(LMIs);

% extract solution variables
Primal.feas = 0;
if p_res > -tol
    Primal.feas = 1;
    Primal.cQ   = cellfun(@double, cQ, 'un', false);
    Primal.cZ   = cellfun(@double, cZ, 'un', false);
    Primal.cL   = cellfun(@double, cL, 'un', false);
    Primal.goal = double(goal);
end

%_______________________________________________________________________
% 2.2: Solve parametric dual

% new LMI system
LMIs = set([]);

% generate the LMI variables
for j=1:get_JNg(A.vertices,p)
    jj = get_monomial_string(KNp(j,:));
    eval(['cU11{' jj '} = sdpvar(nx,nx,''symmetric'');']);
    eval(['cU12{' jj '} = sdpvar(nx,nw,''full'');']);
    eval(['cU22{' jj '} = sdpvar(nw,nw,''symmetric'');']);
    eval(['cU{' jj '} = [cU11{' jj '}, cU12{' jj '}; cU12{' jj '}'', cU22{' jj '}];']);   
    eval(['cV{' jj '} = sdpvar(nz,nz,''symmetric'');']);
end
U11 = HPPDmatrix(cU11);
U12 = HPPDmatrix(cU12);
U22 = HPPDmatrix(cU22);
U   = [U11, U12; U12', U22];
V   = HPPDmatrix(cV);
Iu  = HPPDmatrix(eye(nu));

% construct the LMI terms
W = [A'*U11 + U11*A + Cz'*V*Cz, U11*Bu; Bu'*U11, Iu];

% apply Pólya relaxation of degree d
Ud = U.increase_degree(d);
Vd = V.increase_degree(d);
W  = W.increase_degree(d);

% sufficient LMIs
for j = 1:length(Ud.getcoeffs)
    LMIs = LMIs + set(Ud.getcoeffs{j} > 0);
end
for j = 1:length(Vd.getcoeffs)
    LMIs = LMIs + set(Vd.getcoeffs{j} > 0);
end
for j = 1:length(W.getcoeffs);
    LMIs = LMIs + set(W.getcoeffs{j} > 0);
end

% objective function
goal1 = U12*Bw'; goal2 = U12'*Bw; goal3 = -U22; goal4 = -mu*V;
goal = goal1.trace + goal2.trace + goal3.trace + goal4.trace;
goal = goal.integral;

% solve the SDP
sol = solvesdp(LMIs, -goal, sdpsettings('solver', 'sdpt3', 'sdpt3.maxit', 100, 'dualize', 0, 'verbose', 1));
[p_res,d_res] = checkset(LMIs);

% extract solution variables
Dual.feas = 0;
if p_res > -tol
    Dual.feas = 1;
    Dual.cU   = cellfun(@double, cU, 'un', false);
    Dual.cV   = cellfun(@double, cV, 'un', false);
    Dual.goal = double(goal);
end

%_______________________________________________________________________
% 2.3: Plot results

% sampled solution
figure(1); plot(k2_vec, gamma_vec); hold on;

% primal objective vs parameter
Z = HPPDmatrix(Primal.cZ); 
a1 = linspace(0, 1, 101); a2 = 1-a1;
for j = 1:length(a1) 
    point(j) = Z.trace.f([a1(j) a2(j)]);
end
figure(1); plot(k2L+(k2U-k2L)*a1, point,'r','Linewidth',1.5)

% dual objective vs parameter
goal1 = U12*Bw'; goal2 = U12'*Bw; goal3 = -U22; goal4 = -mu*V;
goal = goal1.trace + goal2.trace + goal3.trace + goal4.trace;
for j = 1:length(a1)  
    point(j) = goal.f([a1(j) a2(j)]);
end
figure(1); plot(k2L+(k2U-k2L)*a1, point,'g','Linewidth',1.5)

% upper/lower bounds of parameter
yL = get(gca,'YLim');
line([k2L k2L],yL,'color','black');
line([k2U k2U],yL,'color','black');

end