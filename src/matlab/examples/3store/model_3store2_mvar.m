% LPV model of 3store example, with parameters (k2,d3)
clear all; close all; clear global; clc

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
     0    , -d3  , d3 ]; % D in terms of parameter d3 
K = [k1+k2, -k2  , 0  ;
     -k2  , k2+k3, -k3;
     0    , -k3  , k3 ]; % K in terms of parameter k2 

E  = blkdiag(eye(3), M);
Ei = blkdiag(eye(3), diag(1./[m1,m2,m3]));
A_nom = Ei*[zeros(3), eye(3); -K, -D]; % nominal value for A matrix
Bu = 1e3 *Ei*[zeros(3); eye(3)];
Bw = Ei*[zeros(3,1); m1; m2; m3];
Cz = 1e3 * [ 1, 0, 0, zeros(1,3);
            -1, 1, 0, zeros(1,3);
             0,-1, 1, zeros(1,3)];

nx = 6; nu = 3; nw = 1; nz = 3;

% lower and upper bounds on parameters (k2,d3)
k2L = k2/2; k2U = 2*k2; d3L = d3/2; d3U = d3;

% define system in terms of Bsplines
kd = parameter(2);  % define (k2,d3) space
kd = kd.to_bspline({[k2L,k2U],[d3L,d3U]});
k2 = kd(1); d3 = kd(2);

D = [d1+d2, -d2  , 0  ;
     -d2  , d2+d3, -d3;
     0    , -d3  , d3 ];           % D in terms of parameter d3 
K = [k1+k2, -k2  , 0  ;
     -k2  , k2+k3, -k3;
     0    , -k3  , k3 ];           % K in terms of parameter k2 
A = Ei*[zeros(3), eye(3); -K, -D]; % A as a function of (k2,d3)

% save('model_3store2_multivar.mat')

