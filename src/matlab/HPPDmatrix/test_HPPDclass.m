% test HPPDmatrix class
clear all
close all
clear global
clc

% 1: create a HPPDmatrix (works!)
coeffs = cell(3,3);
coeffs{3,1} = eye(2); coeffs{2,2} = [1 2; 3 4]; coeffs{1,3} = [1 5; 1 2];
A = HPPDmatrix(coeffs)

% 2: negation of HPPDmatrix (works!)
B = -A;
A.coeffs{1,3}
B.coeffs{1,3}

% 3: add two HPPD matrices 
C = A + B;    % A and B same degree (works!)
A.coeffs{1,3}
B.coeffs{1,3}
C.coeffs{1,3}

B = A*A;     % A and B different degree (works!) 
C = B + A;
A.coeffs{1,end}
B.coeffs{1,end}
C.coeffs{1,end}

% 4: substract two HPPD matrices (works!)
C = A - B;
A.coeffs{1,3}
B.coeffs{1,3}
C.coeffs{1,3}

% 5: take transpose of HPPD matrix (works!)
B = A';
A.coeffs{1,3}
B.coeffs{1,3}

% 6: multiply two HPPD matrices 
B = 2*eye(2);

C = B*B        % two constant matrices (works!)

C = A*B;       % constant matrix and HPPD matrix (works!)
A.coeffs{1,3}
C.coeffs{1,3}

C = A*A;       % two HPPD matrices (works!)
A.coeffs{1,3}
C.coeffs{1,5}

% 7: increase polynomial degree
d = 2;
C = increase_degree(A,d) % HPPD matrix (works!)

% C = increase_degree(B,d,2) % constant matrix


