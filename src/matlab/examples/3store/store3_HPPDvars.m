
clear all
% close all
clear global
clc

%_______________________________________________________________________
% STEP 1. Load model
load('model_3store2.mat');

%_______________________________________________________________________
% STEP 2. Performance ifv k2

% 2.1. Sampled
mu = 15;
k2_vec = linspace(0.5*k2, 2*k2, 100)';
gamma_vec = zeros(size(k2_vec));

% for i = 1:length(k2_vec)
%     k2 = k2_vec(i);
%     K = [k1+k2, -k2  , 0  ;
%          -k2  , k2+k3, -k3;
%          0    , -k3  , k3 ];
%     A = Ei*[zeros(3), eye(3); -K, -D];
%     
%     Q = sdpvar(n);
%     L = sdpvar(nu,n);
%     Z = sdpvar(nu);
%     eps = 1e-9;
% 
%     constr = set(Q >= eps*eye(n)) + set([Q*A'+A*Q'+Bu*L+L'*Bu', Bw; Bw', -eye(nw)] <= 0);
%     constr = constr + set(Cz*Q*Cz' <= mu*eye(nz));
%     constr = constr + set([Q, L'; L, Z] >= 0);
%     sol = solvesdp(constr, trace(Z), sdpsettings('solver', 'sdpt3', 'sdpt3.maxit', 100, 'dualize', 1, 'verbose', 0))
%     gamma_vec(i) = double(trace(Z));
% end
% 
% save gamma_k2 k2_vec gamma_vec
load gamma_k2

%_______________________________________________________________________      
% STEP 2: solve parametric LMI optimization problems

tol = 1e-6; % tolerance primal residual
eps = 0;    % tuning parameter to enforce strict inequality of LMIs

N = vertices(A);                 % #vertices plant
g = degr(A);                     % polynomial degree A(a)
p = 2;                           % degree of Lyapunov matrix Q(a)
d = 0;                           % degree of Pólya relaxation
insknots = []; % linspace(0.1,0.9,20); % inserted knots

% monomial sets and cardinalities
KNg   = get_KNg(N,g);     JNg   = get_JNg(N,g);
KNp   = get_KNg(N,p);     JNp   = get_JNg(N,p);
KNp1  = get_KNg(N,p+1);   JNp1  = get_JNg(N,p+1);
KNpd  = get_KNg(N,p+d);   JNpd  = get_JNg(N,p+d);
KNgp  = get_KNg(N,g+p);   JNgp  = get_JNg(N,g+p);
KNgp1 = get_KNg(N,g+p+1); JNgp1 = get_JNg(N,g+p+1);
KNp1d = get_KNg(N,p+1+d); JNp1d = get_JNg(N,p+1+d);
KNgpd = get_KNg(N,g+p+d); JNgpd = get_JNg(N,g+p+d);

%_______________________________________________________________________ 
% 2.1: Solve parametric primal

% new LMI system
LMIs = set([]);
Primal.rows = 0;

% generate the LMI variables
for j=1:JNp
    jj = gMS(KNp(j,:));
    eval(['Q{' jj '} = sdpvar(nx,nx,''symmetric'');']);
end
for j=1:JNp1
    jj = gMS(KNp1(j,:));
    eval(['Z{' jj '} = sdpvar(nu,nu,''symmetric'');']);
    eval(['L{' jj '} = sdpvar(nu,nx,''full'');']);
end

% construct the LMI terms, all of degree g+p
U11 = HPPD_add(HPPD_tran(HPPD_mult(A,Q)),HPPD_mult(A,Q));
U11 = HPPD_add(U11,HPPD_add(HPPD_mult(Bu,L),HPPD_tran(HPPD_mult(Bu,L))));
U12 = inc_deg(Bw,g+p,N);
U22 = inc_deg(-eye(nw),g+p,N);

V11 = HPPD_mult(HPPD_mult(Cz,Q),HPPD_tran(Cz));
V11 = HPPD_add(V11,inc_deg(-mu*eye(nz),p,N));
V11 = inc_deg(V11,1);

W11 = inc_deg(Q,1);
W12 = HPPD_tran(L);
W22 = Z;

% apply Pólya relaxation of degree d, resulting in terms of degree g+p+d
Qd  = inc_deg(Q,d);
U11 = inc_deg(U11,d);
U12 = inc_deg(U12,d);
U22 = inc_deg(U22,d);
V11 = inc_deg(V11,d);
W11 = inc_deg(W11,d);
W12 = inc_deg(W12,d);
W22 = inc_deg(W22,d);

% determine polynomial coefficients
[coefsU,coefsV,coefsW] = deal(cell(1,JNgpd));
for j = 1:JNgpd
    jj = gMS(KNgpd(j,:));
    coefsU{j} = eval(['[U11{' jj '}  , U12{' jj '};'...
                      ' U12{' jj '}'', U22{' jj '}];']);
    coefsV{j} = eval(['V11{' jj '};']);
    coefsW{j} = eval(['[W11{' jj '}  , W12{' jj '};'...
                      ' W12{' jj '}'', W22{' jj '}];']);
   
    % convert coefficients from HPPD basis to Bspline basis     
    coefsU{j} = coefsU{j}/nchoosek(JNgpd-1,j-1);
    coefsV{j} = coefsV{j}/nchoosek(JNgpd-1,j-1);
    coefsW{j} = coefsW{j}/nchoosek(JNgpd-1,j-1);
end

% knot insertion
degree  = degr(U11);
knots1  = zeros(1,degree+1); 
knots2  = ones(1,degree+1);
knots   = [knots1, knots2];
b       = BSplineBasis(knots,degree); % define BSpline basis
b2      = b.insert_knots(insknots);   % basis after knot insertion
tranmat = b2.transform(b);            % transformation matrix

vecU      = vertcat(coefsU{:});          
vecUnew   = kron(tranmat,eye(nx+nw))*vecU;
newcoefsU = mat2cell(vecUnew,(nx+nw)*ones(1,JNgpd+length(insknots)));

vecV      = vertcat(coefsV{:});          
vecVnew   = kron(tranmat,eye(nz))*vecV;
newcoefsV = mat2cell(vecVnew,(nz)*ones(1,JNgpd+length(insknots)));

vecW      = vertcat(coefsW{:});          
vecWnew   = kron(tranmat,eye(nx+nu))*vecW;
newcoefsW = mat2cell(vecWnew,(nx+nu)*ones(1,JNgpd+length(insknots)));

% sufficient LMIs
for j = 1:JNpd % monomials Lyapunov matrix positive definite
    jj = gMS(KNpd(j,:));
    Term = eval(['Qd{' jj '};']);
    LMIs = LMIs + set(Term > eps*eye(nx));
    Primal.rows = Primal.rows + nx;
end
for j = 1:JNgpd+length(insknots)
    jj = gMS(j-1);
    
    Term = eval(['newcoefsU{' jj '};']);
    LMIs = LMIs + set(Term < -eps*eye(nx+nw));
    Primal.rows = Primal.rows + nx + nw;
    
    Term = eval(['newcoefsV{' jj '};']);
    LMIs = LMIs + set(Term < -eps*eye(nz));
    Primal.rows = Primal.rows + nz;
    
    Term = eval(['newcoefsW{' jj '};']);
    LMIs = LMIs + set(Term > eps*eye(nx+nu));
    Primal.rows = Primal.rows + nx + nu;
end
Primal.vars = size(getvariables(LMIs),2); %number of LMI variables

% objective function
% if the polynomial matrix function Z(a) is expressed in the Bspline basis: 
% int_{0}^{1} trace(Z(a)) da = sum(trace(Z_k))/(degr(Z(a))+1) 
% (use Beta function to derive this property)
coefsZ = cell(1,JNp1);
goal = 0;
for j = 1:JNp1
    jj = gMS(KNp1(j,:));
    coefsZ{j} = eval(['Z{' jj '};']); 
    coefsZ{j} = coefsZ{j}/nchoosek(JNp1-1,j-1); % convert coefs Z(a) from HPPD to Bspline   
    goal = goal + trace(coefsZ{j});
end
goal = goal / (p+2); % degr(Z(a)) = p + 1

% solve the SDP
sol = solvesdp(LMIs, goal, sdpsettings('solver', 'sdpt3', 'sdpt3.maxit', 100, 'dualize', 0, 'verbose', 1))
[p_res,d_res] = checkset(LMIs);

% extract solution variables
Primal.feas = 0;
if p_res > -tol
    Primal.feas = 1;
    for j=1:JNp
        jj = gMS(KNp(j,:));
        eval(['Primal.Q{' jj '} = double(Q{' jj '});'])
    end
    for j=1:JNp1
        jj = gMS(KNp1(j,:));
        eval(['Primal.Z{' jj '} = double(Z{' jj '});'])
        eval(['Primal.L{' jj '} = double(L{' jj '});'])
    end
    Primal.goal = double(goal);
end

% % check the solution
% clear A Q Z L
% a1_vec = linspace(0,1,1e2)';
% check_LMI1 = zeros(size(a1_vec));
% check_LMI2 = zeros(size(a1_vec));
% check_LMI3 = zeros(size(a1_vec));
% for i = 1:length(a1_vec)
%     a1 = a1_vec(i);
%     a2 = 1-a1;
%     k2 = k2U*a1+k2L*a2;
%     K = [k1+k2, -k2  , 0  ;
%         -k2  , k2+k3, -k3;
%         0    , -k3  , k3 ];
%     A = Ei*[zeros(3), eye(3); -K, -D];
%     Q = Info.Q{3,1}*a1^2 + Info.Q{2,2}*a1*a2 + Info.Q{1,3}*a2^2;
%     Z = Info.Z{4,1}*a1^3 + Info.Z{3,2}*a1^2*a2 + Info.Z{2,3}*a1*a2^2 + Info.Z{1,4}*a2^3;
%     L = Info.L{4,1}*a1^3 + Info.L{3,2}*a1^2*a2 + Info.L{2,3}*a1*a2^2 + Info.L{1,4}*a2^3;
%     
%     LMI1 = [Q*A'+A*Q+Bu*L+L'*Bu', Bw; Bw', -eye(nw)];
%     LMI2 = Cz*Q*Cz'-mu*eye(nz);
%     LMI3 = [Q, L'; L, Z];
%     
%     check_LMI1(i) = max(eig(LMI1));
%     check_LMI2(i) = max(eig(LMI2));
%     check_LMI3(i) = min(eig(LMI3));
% end

% figure
% subplot(131)
% plot(check_LMI1)
% subplot(132)
% plot(check_LMI2)
% subplot(133)
% plot(check_LMI3)

% objective function as a function of the parameter
a1_vec = (k2_vec-k2L) / (k2U-k2L);
a2_vec = (k2U-k2_vec) / (k2U-k2L);
gamma_vec2 = 0;
for j=1:JNp1
    gamma_vec2 = gamma_vec2 + trace(Primal.Z{JNp1-j+1,j})*a1_vec.^(JNp1-j).*a2_vec.^(j-1);
end

%_______________________________________________________________________
% 2.2: Solve parametric dual
clear U11 U12 U22
 
% new LMI system
LMIs = set([]);
Dual.rows = 0;

% generate the LMI variables
for j=1:JNp
    jj = gMS(KNp(j,:));
    eval(['U11{' jj '} = sdpvar(nx,nx,''symmetric'');']);
    eval(['U12{' jj '} = sdpvar(nx,nw,''full'');']);
    eval(['U22{' jj '} = sdpvar(nw,nw,''symmetric'');']);
    eval(['V{' jj '} = sdpvar(nz,nz,''symmetric'');']);
end

% construct the LMI terms
T11 = HPPD_add(HPPD_tran(HPPD_mult(U11,A)),HPPD_mult(U11,A));
T11 = HPPD_add(T11,inc_deg(HPPD_mult(HPPD_mult(HPPD_tran(Cz),V),Cz),1));
T12 = inc_deg(HPPD_mult(U11,Bu),1);
T22 = inc_deg(eye(nu),g+p,N);

% apply Pólya relaxation of degree d
U11d = inc_deg(U11,d);
U12d = inc_deg(U12,d);
U22d = inc_deg(U22,d);
Vd   = inc_deg(V,d);
T11  = inc_deg(T11,d);
T12  = inc_deg(T12,d);
T22  = inc_deg(T22,d);

% determine polynomial coefficients
[coefsU,coefsV] = deal(cell(1,JNpd));
for j = 1:JNpd
    jj = gMS(KNpd(j,:));
    coefsU{j} = eval(['[U11d{' jj '}  , U12d{' jj '};'...
                      ' U12d{' jj '}'', U22d{' jj '}];']);
    coefsV{j} = eval(['Vd{' jj '};']);
   
    % convert coefficients from HPPD basis to Bspline basis     
    coefsU{j} = coefsU{j}/nchoosek(JNpd-1,j-1);
    coefsV{j} = coefsV{j}/nchoosek(JNpd-1,j-1);
end
coefsT = cell(1,JNgpd);
for j = 1:JNgpd
    jj = gMS(KNgpd(j,:));
    coefsT{j} = eval(['[T11{' jj '}  , T12{' jj '};'...
                      ' T12{' jj '}'', T22{' jj '}];']);
   
    % convert coefficients from HPPD basis to Bspline basis     
    coefsT{j} = coefsT{j}/nchoosek(JNgpd-1,j-1);
end

% knot insertion
degree  = degr(U11d);
knots1  = zeros(1,degree+1); 
knots2  = ones(1,degree+1);
knots   = [knots1, knots2];
b       = BSplineBasis(knots,degree); % define BSpline basis
b2      = b.insert_knots(insknots);   % basis after knot insertion
tranmat = b2.transform(b);            % transformation matrix

vecU    = vertcat(coefsU{:});          
vecUnew = kron(tranmat,eye(nx+nw))*vecU;
newU    = mat2cell(vecUnew,(nx+nw)*ones(1,JNpd+length(insknots)));

vecV    = vertcat(coefsV{:});          
vecVnew = kron(tranmat,eye(nz))*vecV;
newV    = mat2cell(vecVnew,nz*ones(1,JNpd+length(insknots)));

degree  = degr(T11);
knots1  = zeros(1,degree+1); 
knots2  = ones(1,degree+1);
knots   = [knots1, knots2];
b       = BSplineBasis(knots,degree); % define BSpline basis
b2      = b.insert_knots(insknots);   % basis after knot insertion
tranmat = b2.transform(b);            % transformation matrix

vecT    = vertcat(coefsT{:});          
vecTnew = kron(tranmat,eye(nx+nu))*vecT;
newT    = mat2cell(vecTnew,(nx+nu)*ones(1,JNgpd+length(insknots)));

% sufficient LMIs
for j = 1:JNpd+length(insknots)
    jj = gMS(j-1);
    Term = eval(['newU{' jj '};']);
    LMIs = LMIs + set(Term > eps*eye(nx+nw));
    Dual.rows = Dual.rows + nx + nw;
    Term = eval(['newV{' jj '};']);
    LMIs = LMIs + set(Term > eps*eye(nz));
    Dual.rows = Dual.rows + nz;
end
for j = 1:JNgpd+length(insknots)
    jj = gMS(j-1);
    Term = eval(['newT{' jj '};']);
    LMIs = LMIs + set(Term > eps*eye(nx+nu));
    Dual.rows = Dual.rows + nx + nu;
end
Dual.vars = size(getvariables(LMIs),2); %number of LMI variables

% objective function
goal = 0;
for j = 1:JNpd
    jj = gMS(j-1);
    goal = goal + eval(['trace(coefsU{' jj '}(1:nx,nx+1:end)*Bw'')+trace(coefsU{' jj '}(1:nx,nx+1:end)''*Bw)-trace(coefsU{' jj '}(nx+1:end,nx+1:end))-trace(mu*coefsV{' jj '});']);
end
goal = -goal / (p+d+1); % degr(Ud(a)) = degr(Vd(a)) = p + d

% solve the SDP
sol = solvesdp(LMIs, goal, sdpsettings('solver', 'sdpt3', 'sdpt3.maxit', 100, 'dualize', 0, 'verbose', 1))
[p_res,d_res] = checkset(LMIs);

% extract solution variables
Dual.feas = 0;
if p_res > -tol
    Dual.feas = 1;
    for j=1:JNp
        jj = gMS(KNp(j,:));
        eval(['Dual.U11{' jj '} = double(U11{' jj '});'])
        eval(['Dual.U12{' jj '} = double(U12{' jj '});'])
        eval(['Dual.U22{' jj '} = double(U22{' jj '});'])
        eval(['Dual.V{' jj '} = double(V{' jj '});'])
    end
    Dual.goal = -double(goal);
end

% % % check the solution
% clear A 
% a1_vec = linspace(0,1,1e2)';
% check_LMI1 = zeros(size(a1_vec));
% check_LMI2 = zeros(size(a1_vec));
% check_LMI3 = zeros(size(a1_vec));
% for i = 1:length(a1_vec)
%     a1 = a1_vec(i);
%     a2 = 1-a1;
%     k2 = k2U*a1+k2L*a2;
%     K = [k1+k2, -k2  , 0  ;
%         -k2  , k2+k3, -k3;
%         0    , -k3  , k3 ];
%     A = Ei*[zeros(3), eye(3); -K, -D];
%     U11 = Info.U11{3,1}*a1^2 + Info.U11{2,2}*a1*a2 + Info.U11{1,3}*a2^2;
%     U12 = Info.U12{3,1}*a1^2 + Info.U12{2,2}*a1*a2 + Info.U12{1,3}*a2^2;
%     U22 = Info.U22{3,1}*a1^2 + Info.U22{2,2}*a1*a2 + Info.U22{1,3}*a2^2;
%     V = Info.V{3,1}*a1^2 + Info.V{2,2}*a1*a2 + Info.V{1,3}*a2^2;
%      
%     LMI1 = [A'*U11+U11*A+Cz'*V*Cz, U11*Bu; Bu'*U11, eye(nu)];
%     LMI2 = [U11, U12; U12', U22];
%     LMI3 = V;
%     
%     check_LMI1(i) = min(eig(LMI1));
%     check_LMI2(i) = min(eig(LMI2));
%     check_LMI3(i) = min(eig(LMI3));
% end
% 
% % figure
% % subplot(131)
% % plot(check_LMI1)
% % subplot(132)
% % plot(check_LMI2)
% % subplot(133)
% % plot(check_LMI3)

% objective function as a function of the parameter
gamma_vec3 = 0;
for j=1:JNp
    gamma_vec3 = gamma_vec3 + a1_vec.^(JNp-j).*a2_vec.^(j-1)*(...
        trace(Dual.U12{JNp-j+1,j}*Bw') + ...
        trace(Dual.U12{JNp-j+1,j}'*Bw) - ...
        trace(Dual.U22{JNp-j+1,j})     - ...
        trace(mu*Dual.V{JNp-j+1,j}));
end

% figure
plot(k2_vec, gamma_vec,'LineWidth',1.5); 
hold on;
plot(k2_vec, gamma_vec2,'LineWidth',1.5,'color','red');
plot(k2_vec, gamma_vec3,'g--','LineWidth',1.5);

yL = get(gca,'YLim');
line([k2L k2L],yL,'color','black');
line([k2U k2U],yL,'color','black');