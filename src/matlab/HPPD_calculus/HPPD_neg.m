function C = HPPD_neg(A) 
% Negation of a HPPD matrix A(a), 
% resulting in C(a) = -A(a) 

N   = vertices(A);  % #vertices of A
g   = degr(A);      % polynomial degree of HPPD matrix
KNg = get_KNg(N,g); % monomial set
JNg = get_JNg(N,g); % cardinality 

for j = 1:JNg
    jj = KNg(j,:);
    eval(['C{' gMS(jj) '} = -A{' gMS(jj) '};'])
end