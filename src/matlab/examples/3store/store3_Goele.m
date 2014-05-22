
clear all
close all
clear global
clc



% STEP 1. System data
% -------------------

m1 = 5.897e3;       % [kg]
m2 = 5.897e3;       % [kg]
m3 = 5.897e3;       % [kg]
k1 = 33.732e6;      % [N/m]
k2 = 29.093e6;      % [N/m]
k3 = 28.621e6;      % [N/m]
d1 = 67e3;          % [Ns/m]
d2 = 58e3;          % [Ns/m]
d3 = 57e3;          % [Ns/m]

d2 = 2*d2;          % -> more interesting for opt. k2

M = diag([m1,m2,m3]);
D = [d1+d2, -d2  , 0  ;
     -d2  , d2+d3, -d3;
     0    , -d3  , d3 ];
K = [k1+k2, -k2  , 0  ;
     -k2  , k2+k3, -k3;
     0    , -k3  , k3 ];
 
E = blkdiag(eye(3), M);
Ei = blkdiag(eye(3), diag(1./[m1,m2,m3]));
A = Ei*[zeros(3), eye(3); -K, -D];
Bu = 1e3*Ei*[zeros(3); eye(3)];
Bw = Ei*[zeros(3,1); m1; m2; m3];
Cz = 1e3 * [ 1, 0, 0, zeros(1,3);
            -1, 1, 0, zeros(1,3);
             0,-1, 1, zeros(1,3)];

n = 6;
nu = 3;
nw = 1;
nz = 3;



% STEP 2. Performance ifv k2
% --------------------------

% 2.1. Sampled
mu = 15;
k2_vec = linspace(0.5*k2, 2*k2, 100)';
% gamma_vec = zeros(size(k2_vec));

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

% save gamma_k2 k2_vec gamma_vec
load gamma_k2


% 2.2. Parametric primal
K0 = [k1, 0  , 0  ;
      0 , k3 , -k3;
      0 , -k3, k3 ];
K1 = [1, -1, 0; -1, 1, 0; 0, 0, 0];  
A0 = Ei*[zeros(3), eye(3); -K0, -D];
A1 = Ei*[zeros(3), zeros(3); -K1, zeros(3)];

k2L = 2.5e7;%k2;
k2U = 3.5e7;%  2*k2;
eps = 0;%1e-7;

Q20 = sdpvar(n);
Q11 = sdpvar(n);
Q02 = sdpvar(n);
Z30 = sdpvar(nu);
Z21 = sdpvar(nu);
Z12 = sdpvar(nu);
Z03 = sdpvar(nu);
L30 = sdpvar(nu,n);
L21 = sdpvar(nu,n);
L12 = sdpvar(nu,n);
L03 = sdpvar(nu,n);

U30 = (Q20*A0'+A0*Q20') + k2L*(Q20*A1'+A1*Q20') + Bu*L30 + L30'*Bu';
U21 = (Q20*A0'+A0*Q20') + (Q11*A0'+A0*Q11') + Bu*L21 + L21'*Bu'...
       + k2L*(Q11*A1'+A1*Q11') + k2U*(Q20*A1'+A1*Q20');
U12 = (Q02*A0'+A0*Q02') + (Q11*A0'+A0*Q11') + Bu*L12 + L12'*Bu'...
       + k2L*(Q02*A1'+A1*Q02') + k2U*(Q11*A1'+A1*Q11');
U03 = (Q02*A0'+A0*Q02') + k2U*(Q02*A1'+A1*Q02') + Bu*L03 + L03'*Bu';
   
constr = set(Q20 >= eps*eye(n)) + set(Q11 >= eps*eye(n)) + set(Q02 >= eps*eye(n));
constr = constr + set([U30, Bw; Bw', -eye(nw)] <= -eps) ... 
                + set([U21, 3*Bw; 3*Bw', -3*eye(nw)] <= -eps) ... 
                + set([U12, 3*Bw; 3*Bw', -3*eye(nw)] <= -eps) ... 
                + set([U03, Bw; Bw', -eye(nw)] <= -eps) ;
constr = constr + set(Cz*Q20*Cz' <= mu*eye(nz)) ...
                + set(Cz*(Q20+Q11)*Cz' <= 3*mu*eye(nz)) ...
                + set(Cz*(Q11+Q02)*Cz' <= 3*mu*eye(nz)) ...
                + set(Cz*Q02*Cz' <= mu*eye(nz)) ;            
constr = constr + set([Q20, L30'; L30, Z30] >= eps) ...
                + set([Q20+Q11, L21'; L21, Z21] >= eps) ...
                + set([Q11+Q02, L12'; L12, Z12] >= eps) ...
                + set([Q02, L03'; L03, Z03] >= eps);         
goal = 1/4*trace(Z30) + 1/12*trace(Z21) + 1/12*trace(Z12) + 1/4*trace(Z03);            

sol = solvesdp(constr, goal, sdpsettings('solver', 'sdpt3', 'sdpt3.maxit', 100, 'dualize', 0, 'verbose', 1))
Q20 = double(Q20);   Q11 = double(Q11);   Q02 = double(Q02);
Z30 = double(Z30);   Z21 = double(Z21);   Z12 = double(Z12);   Z03 = double(Z03);
L30 = double(L30);   L21 = double(L21);   L12 = double(L12);   L03 = double(L03);

l1_vec = linspace(0,1,1e2)';
check_LMI1 = zeros(size(l1_vec));
check_LMI2 = zeros(size(l1_vec));
check_LMI3 = zeros(size(l1_vec));
for i = 1:length(l1_vec)
   l1 = l1_vec(i);
   l2 = 1-l1;
   k2 = k2L*l1+k2U*l2;
   K = [k1+k2, -k2  , 0  ;
        -k2  , k2+k3, -k3;
        0    , -k3  , k3 ];
   A = Ei*[zeros(3), eye(3); -K, -D];
   Q = Q20*l1^2 + Q11*l1*l2 + Q02*l2^2;
   Z = Z30*l1^3 + Z21*l1^2*l2 + Z12*l1*l2^2 + Z03*l2^3;
   L = L30*l1^3 + L21*l1^2*l2 + L12*l1*l2^2 + L03*l2^3;
   
   LMI1 = [Q*A'+A*Q+Bu*L+L'*Bu', Bw; Bw', -eye(nw)];
   LMI2 = Cz*Q*Cz'-mu*eye(nz);
   LMI3 = [Q, L'; L, Z];
   
   check_LMI1(i) = max(eig(LMI1));
   check_LMI2(i) = max(eig(LMI2));
   check_LMI3(i) = min(eig(LMI3));
end

figure
subplot(131)
plot(check_LMI1)
subplot(132)
plot(check_LMI2)
subplot(133)
plot(check_LMI3)


l1_vec = (k2U-k2_vec) / (k2U-k2L);
l2_vec = (k2_vec-k2L) / (k2U-k2L);
gamma_vec2 = trace(Z30)*l1_vec.^3 + trace(Z21)*l1_vec.^2.*l2_vec + trace(Z12)*l1_vec.*l2_vec.^2 + trace(Z03)*l2_vec.^3;

figure
plot(k2_vec, [gamma_vec, gamma_vec2])


% 2.3. Parametric dual
U11_20 = sdpvar(n);
U11_11 = sdpvar(n);
U11_02 = sdpvar(n);
U12_20 = sdpvar(n, nw, 'full');
U12_11 = sdpvar(n, nw, 'full');
U12_02 = sdpvar(n, nw, 'full');
U22_20 = sdpvar(nw);
U22_11 = sdpvar(nw);
U22_02 = sdpvar(nw);
V_20 = sdpvar(nz);
V_11 = sdpvar(nz);
V_02 = sdpvar(nz);

Q30 = (A0'*U11_20+U11_20*A0) + Cz'*V_20*Cz ...
        + k2L*(A1'*U11_20+U11_20*A1);
Q21 = (A0'*U11_20+U11_20*A0) + (A0'*U11_11+U11_11*A0) + Cz'*V_20*Cz + Cz'*V_11*Cz...
        + k2L*(A1'*U11_11+U11_11*A1) + k2U*(A1'*U11_20+U11_20*A1);
Q12 = (A0'*U11_11+U11_11*A0) + (A0'*U11_02 + U11_02*A0) + Cz'*V_11*Cz + Cz'*V_02*Cz...
        + k2L*(A1'*U11_02+U11_02*A1) + k2U*(A1'*U11_11+U11_11*A1);
Q03 = (A0'*U11_02 + U11_02*A0) + Cz'*V_02*Cz  ...
        + k2U*(A1'*U11_02+U11_02*A1);

constr = set([U11_20, U12_20; U12_20', U22_20] >= eps*eye(n+nw)) ...
       + set([U11_11, U12_11; U12_11', U22_11] >= eps*eye(n+nw)) ...
       + set([U11_02, U12_02; U12_02', U22_02] >= eps*eye(n+nw)) ;
constr = constr + set(V_20 >= eps*eye(nz)) ...
                + set(V_11 >= eps*eye(nz)) ...
                + set(V_02 >= eps*eye(nz)) ;
constr = constr + set([Q30, U11_20*Bu; Bu'*U11_20, eye(nu)] >= 0) ...
                + set([Q21, (U11_11+U11_20)*Bu; Bu'*(U11_11+U11_20), 3*eye(nu)] >= 0) ...
                + set([Q12, (U11_02+U11_11)*Bu; Bu'*(U11_02+U11_11), 3*eye(nu)] >= 0) ...
                + set([Q03, (U11_02)*Bu; Bu'*(U11_02), eye(nu)] >= 0) ;
goal = - 1/3* (trace(U12_20*Bw') + trace(U12_20'*Bw) - trace(U22_20) - trace(mu*V_20)) ...
       - 1/6* (trace(U12_11*Bw') + trace(U12_11'*Bw) - trace(U22_11) - trace(mu*V_11)) ...
       - 1/3* (trace(U12_02*Bw') + trace(U12_02'*Bw) - trace(U22_02) - trace(mu*V_02)) ;

sol = solvesdp(constr, goal, sdpsettings('solver', 'sdpt3', 'sdpt3.maxit', 100, 'dualize', 0, 'verbose', 1))

U11_20 = double(U11_20);   U11_11 = double(U11_11);   U11_02 = double(U11_02);
U12_20 = double(U12_20);   U12_11 = double(U12_11);   U12_02 = double(U12_02);
U22_20 = double(U22_20);   U22_11 = double(U22_11);   U22_02 = double(U22_02);
V_20 = double(V_20);       V_11 = double(V_11);       V_02 = double(V_02);

gamma_vec3 =  l1_vec.^2      * (trace(U12_20*Bw') + trace(U12_20'*Bw) - trace(U22_20) - trace(mu*V_20)) + ...
              l1_vec.*l2_vec * (trace(U12_11*Bw') + trace(U12_11'*Bw) - trace(U22_11) - trace(mu*V_11)) + ...
              l2_vec.^2      * (trace(U12_02*Bw') + trace(U12_02'*Bw) - trace(U22_02) - trace(mu*V_02)) ;
          
figure
plot(k2_vec, [gamma_vec, gamma_vec2, gamma_vec3])          