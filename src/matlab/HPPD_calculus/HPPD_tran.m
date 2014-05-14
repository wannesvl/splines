function C = HPPD_tran(A)
% determine the transpose of a HPPD matrix A(a),
% resulting in C(a) = A(a)^T

if ~iscell(A) % if A is not a cell (but a constant matrix)
    C = A';
else
    N   = vertices(A);  % #vertices A(a)
    g   = degr(A);      % polynomial degree A(a)
    KNg = get_KNg(N,g); % monomial set
    JNg = get_JNg(N,g); % cardinality
    
    for j = 1:JNg
        jj = KNg(j,:);
        eval(['C{' gMS(jj) '} = A{' gMS(jj) '}'';'])
    end    
end

end