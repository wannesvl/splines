function C = HPPD_mult(A,B)
% Multiply two HPPD matrices A(a) and B(a), 
% resulting in C(a) = A(a)B(a) 
% 
% The cases that A and/or B are constant (i.e. parameter-independent) 
% matrices are taken into account

if ~iscell(A) && ~iscell(B)     % A and B constant matrices
    type = 1;
elseif ~iscell(A) && iscell(B)  % A constant, B HPPD
    type = 2;
elseif iscell(A) && ~iscell(B)  % A HPPD, B constant
    type = 3;
else                            % A HPPD, B HPPD
    type = 4;
end

switch type

    case 1 % A and B constant matrices
        C = A*B;
        
    case 2 % A constant, B HPPD
        N = vertices(B);
        g = degr(B);
        
        ncA = size(A,2);        
        ind = find(~cellfun('isempty',B),1);   
        nrB = size(B{ind},1); % #columns of B(a)
        if ncA ~= nrB
            error('dimensions of matrices must agree')
        end
        
        % monomial sets and cardinalities
        KNg = get_KNg(N,g); JNg = get_JNg(N,g);
       
        % construct monomials of A*B(a)
        for j = 1:JNg
            jj = gMS(KNg(j,:));
            eval(['C{' jj '} = A*B{' jj '};'])
        end
        
    case 3 % A HPPD, B constant
        N = vertices(A);
        g = degr(A);
     
        ind = find(~cellfun('isempty',A),1);   
        ncA = size(A{ind},2); % #columns of B(a)
        nrB = size(B,1);
        if ncA ~= nrB
            error('dimensions of matrices must agree')
        end
        
        % monomial sets and cardinalities
        KNg = get_KNg(N,g); JNg = get_JNg(N,g);
       
        % construct monomials of A*B(a)
        for j = 1:JNg
            jj = gMS(KNg(j,:));
            eval(['C{' jj '} = A{' jj '}*B;'])
        end     
        
    case 4 % A HPPD, B HPPD

        % dimensions
        N = vertices(A); % #vertices A(a)
        if vertices(B) ~= N
            error('A and B should have same number of vertices!')
        end
        gA  = degr(A); % degree of A(a)
        gB  = degr(B); % degree of B(a)

        ind = find(~cellfun('isempty',A),1);   
        nrA = size(A{ind},1); % #rows of A(a)
        ncA = size(A{ind},2); % #columns of A(a)

        ind = find(~cellfun('isempty',B),1);   
        nrB = size(B{ind},1); % #rows of B(a)
        ncB = size(B{ind},2); % #columns of B(a)
        if ncA ~= nrB
            error('dimensions of matrices must agree')
        end

        % monomial sets and cardinalities
        KNc = get_KNg(N,gA+gB); JNc = get_JNg(N,gA+gB);
        KNb = get_KNg(N,gB);    JNb = get_JNg(N,gB);

        % construct monomials of A(a)B(a)
        Cmon = zeros(nrA,ncB);
        for i=1:JNc
            j = KNc(i,:);
            for m=1:JNb
                k = KNb(m,:);
                if compk(j,k) % if j>=k elementwise    
                    Cmon = Cmon + eval(['A{' gMS(j-k) '}*B{' gMS(k) '};']);
                end
            end
            eval(['C{' gMS(j) '} = Cmon;']);
            Cmon = zeros(nrA,ncB);
        end
end

end
