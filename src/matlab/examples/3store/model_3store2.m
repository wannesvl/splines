% HPPD model of 3store example

clear all
close all
clear global
clc

% system parameters
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
% A = Ei*[zeros(3), eye(3); -K, -D];
Bu = 1e3*Ei*[zeros(3); eye(3)];
Bw = Ei*[zeros(3,1); m1; m2; m3];
Cz = 1e3 * [ 1, 0, 0, zeros(1,3);
            -1, 1, 0, zeros(1,3);
             0,-1, 1, zeros(1,3)];

nx = 6;
nu = 3;
nw = 1;
nz = 3;

% HPPD LPV model with parameter k2
% a1 = (k2-k2L)/(k2U-k2L), a2 = (k2U - k2)/(k2U-k2L)
K0 = [k1, 0  , 0  ;
      0 , k3 , -k3;
      0 , -k3, k3 ];
K1 = [1, -1, 0; -1, 1, 0; 0, 0, 0];  
A0 = Ei*[zeros(3), eye(3); -K0, -D];
A1 = Ei*[zeros(3), zeros(3); -K1, zeros(3)];

k2L = 1.5e7;% k2;
k2U = 5.8e7;% 2*k2;

% A(a1,a2) = a1 A{2,1} + a2 A{1,2}
A{2,1} = A0 + k2U*A1;
A{1,2} = A0 + k2L*A1;

save('model_3store2.mat')

