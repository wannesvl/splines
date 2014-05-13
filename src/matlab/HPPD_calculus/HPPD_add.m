function C = HPPD_add(A,B)
% Add two HPPD matrices A(a) and B(a), 
% resulting in C(a) = A(a) + B(a) 

% dimensions
N = vertices(A);  % #vertices A(a)
if vertices(B) ~= N
    error('A and B should have same number of vertices!')
end

gA = degr(A); % degree of A(a)
gB = degr(B); % degree of B(a)
if gA ~= gB
    error('A and B should have the same polynomial degree')
end

ind = find(~cellfun('isempty',A),1);   
nrA = size(A{ind},1); % #rows of A(a)
ncA = size(A{ind},2); % #columns of A(a)

ind = find(~cellfun('isempty',B),1); 
nrB = size(B{ind},1); % #rows of B(a)
ncB = size(B{ind},2); % #columns of B(a)

if nrA ~= nrB || ncA ~= ncB
    error('dimensions of matrices must agree')
end

% monomial sets and cardinalities
KNg = get_KNg(N,gA); JNg = get_JNg(N,gA);

% construct monomials of A(a) + B(a)
for i=1:JNg
    j = KNg(i,:);
    eval([ 'C{' gMS(j) '} = A{' gMS(j) '} + B{' gMS(j) '};']);
end

end