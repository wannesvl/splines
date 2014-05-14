function Ad = inc_deg(A,d,N)
% Increase the degree of a HPPD matrix A(a) 
% The degree of an HPPD matrix A(a) is increased by d
% 
% The number of vertices N does not need to be specified if A(a) has a
% polynomial degree N >= 2.

if d == 0 % if no increase in degree
    Ad = A;
    return
end

if ~iscell(A) % if A is not a cell (but a constant matrix)
    if nargin < 3
        error('number of vertices unknown');
    end
    
    % monomial sets
    KNd = get_KNg(N,d); JNd = get_JNg(N,d);
    
    % increase degree of constant matrix
    for i=1:JNd
        j  = KNd(i,:); 
        Cd = factorial(d)/Pi(j);     
        eval(['Ad{' gMS(j) '} = Cd*A;'])
    end
else
    N = vertices(A); % #vertices A(a)
    g = degr(A);     % polynomial degree A(a)
    
    ind = find(~cellfun('isempty',A),1);
    nr = size(A{ind},1); % #rows of A(a)
    nc = size(A{ind},2); % #columns of A(a)
    
    % monomial sets and cardinalities
    KNgd = get_KNg(N,g+d); JNgd = get_JNg(N,g+d);
    KNd  = get_KNg(N,d);   JNd  = get_JNg(N,d);

    % increase degree of A(a) to g+d  
    for i=1:JNgd
        Amon = zeros(nr,nc);
        j = KNgd(i,:);    
        for m=1:JNd
            k = KNd(m,:);        
            if compk(j,k) % if j>=k elementwise
                Cd   = factorial(d)/Pi(k);
                Amon = Amon + Cd*eval(['A{' gMS(j-k) '};']);
            end
        end
        eval(['Ad{' gMS(j) '} = Amon;'])
    end
end

end