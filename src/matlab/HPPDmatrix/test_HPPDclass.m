% test HPPDmatrix class
clear all
close all
clear global
clc

% 1: create a HPPDmatrix (works!)
coeffs = cell(3,3);
coeffs{3,1} = eye(2); coeffs{2,2} = [1 2; 3 4]; coeffs{1,3} = [1 5; 1 2];
A = HPPDmatrix(coeffs);

coeffs = [1 2; 0 1]; % constant matrix
B = HPPDmatrix(coeffs);

% 2: negation of HPPDmatrix (works!)
B = -A;

% 3: add two HPPD matrices 
C = A + B;   % A and B same degree (works!)

B = A*A;     % A and B different degree (works!) 
C = B + A;

% 4: substract two HPPD matrices (works!)
C = A - B;

% 5: take transpose of HPPD matrix (works!)
B = A';

% 6: multiply two HPPD matrices 
C = B*B; % two constant matrices (in HPPD class) (works!)
C = A*B; % constant matrix and HPPD matrix (works!)
C = A*A; % two HPPD matrices (works!)

% 7: increase polynomial degree
d = 2;
C = A.increase_degree(d);   % HPPD matrix (works!)
C = B.increase_degree(d,2); % constant matrix (works!)

% 8: concatenation
C = [A, B]; % horizontal concatenation (works!)
C = [A; B]; % vertical concatenation (works!)

O22 = HPPDmatrix(zeros(2,2)); % combined (works!)
I2 = HPPDmatrix(eye(2));
C = [A, O22; I2, B];

% 9: get coefficients (works!)
coeffs = A.getcoeffs;

% 10: evaluation of HPPD matrix at point in unit simplex (works!)
s = A.f([1/2 1/2]);

% 11: test with sdp variables (works!)
c{1,2} = sdpvar(2,2);
c{2,1} = sdpvar(2,2);
A = HPPDmatrix(c);

% 12: trace of HPPD matrix
A.trace;

% 13: integral of HPPD matrix
A.integral;




