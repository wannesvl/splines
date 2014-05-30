classdef HPPDmatrix
    properties (Access=protected)
        cl
    end
    properties
        vertices
        degree
        coeffs
    end
    methods
        function A = HPPDmatrix(coeffs)
            % Constructor for HPPDmatrix
            %
            % Assumptions: - #vertices >= 2
            %              - parameter in unit simplex
            %
            % coeffs is a cell array where the element
            % coeffs{k_1+1,...,k_N+1} corresponds to monomial 
            % k = (k_1,...,k_N), where k_1 + ... + k_N = vertices.
            % All other cell elements are assumed to be empty.
            %
            % vertices and degree are derived from coeffs.
            %
            % Ex.: A(a1,a2) = a1 A_10 + a2 A_01
            % is stored in a 2x2 cell as A{2,1} = A10, A{1,2} = A01,
            % and represents a HPPD matrix of degree 1 with 2 vertices
            A.vertices = length(size(coeffs));
            A.degree   = size(coeffs,1) - 1;
            % TO DO: build in check for equal size coefficients
            A.coeffs   = coeffs;
            A.cl       = str2func(class(A));
        end
 
        function KNg = get_KNg(N,g) 
        %generate all the solutions of: k1 + k2 + ... + kN = g
            KNg = [];
            if N == 1
                KNg = g;
            else
                for i = 0:g
                    temp = get_KNg(N-1,g-i);
                    KNg = [ KNg, i*ones(size(temp,1),1), temp ];
                end
            end
        end
        
        function JNg = get_JNg(N,g) % cardinality of KNg
            JNg = factorial(N+g-1)/(factorial(g)*factorial(N-1));
        end
        
        function b = compk(k1,k2) % check if k1 >= k2 elementwise
            if length(k1)==length(k2)
                b=0;
                if sum(k1>=k2)==length(k1),
                    b=1;
                end
            else
                error('Cannot compare monomial coefficients of different length!')
            end
        end
        
        function mStr = gMS(k) % get monomial string from vector
            k = k+1;
            mStr = num2str(k(1));
            if length(k)>1
                for i = 2:length(k)
                    mStr = [ mStr ',' num2str(k(i)) ];
                end
            end
        end
        
        function R = Pi(monomial) % Pi-coefficient
            R=prod(factorial(monomial));
        end
        
        function s = plus(A,B) 
            if A.vertices ~= B.vertices
                error('self and other should have same number of vertices!')
            end
            if A.degree ~= B.degree % homogenize s.t. A and B have same degree
                if A.degree > B.degree
                   B = increase_degree(B,A.degree-B.degree); 
                else
                   A = increase_degree(A,B.degree-A.degree);  
                end
            end
            ind = find(~cellfun('isempty',A.coeffs),1);   
            nrA = size(A.coeffs{ind},1); 
            ncA = size(A.coeffs{ind},2);
            ind = find(~cellfun('isempty',B.coeffs),1); 
            nrB = size(B.coeffs{ind},1);
            ncB = size(B.coeffs{ind},2); 
            if nrA ~= nrB || ncA ~= ncB
                error('dimensions of matrices must agree')
            end
            coeffs = cellfun(@plus, A.coeffs, B.coeffs, 'UniformOutput', false);
            s = A.cl(coeffs);
        end
        
        function C = uminus(A)
            coeffs = cellfun(@uminus, A.coeffs, 'UniformOutput', false);
            C = A.cl(coeffs);
        end      
        
        function C = minus(A,B)
            C = A + (-B);
        end
        
        function B = ctranspose(A)
            coeffs = cellfun(@ctranspose, A.coeffs, 'UniformOutput', false);
            B = A.cl(coeffs);
        end      
        
        function C = mtimes(A,B)       
            % The case that A or B is a constant (parameter-independent) 
            % matrix is taken into account        
            if isa(A,'double') && isa(B.coeffs,'cell')
                ind = find(~cellfun('isempty',B.coeffs),1);   
                if size(A,2) ~= size(B.coeffs{ind},1)
                    error('dimensions of matrices must agree')
                end 
                KNg = get_KNg(B.vertices,B.degree);                
                for j = 1:get_JNg(B.vertices,B.degree)
                    jj = gMS(KNg(j,:));
                    eval(['coeffs{' jj '} = A*B.coeffs{' jj '};'])
                end
                C = B.cl(coeffs); 
            elseif isa(A.coeffs,'cell') && isa(B,'double')
                ind = find(~cellfun('isempty',A.coeffs),1);  
                if size(A.coeffs{ind},2) ~= size(B,1)
                    error('dimensions of matrices must agree')
                end                   
                KNg = get_KNg(A.vertices,A.degree);                
                for j = 1:get_JNg(A.vertices,A.degree)
                    jj = gMS(KNg(j,:));
                    eval(['coeffs{' jj '} = A.coeffs{' jj '}*B;'])
                end
                C = A.cl(coeffs);  
            else
                if A.vertices ~= B.vertices
                    error('A and B should have same number of vertices!')
                end
                ind = find(~cellfun('isempty',A.coeffs),1);
                nrA = size(A.coeffs{ind},1); 
                ncA = size(A.coeffs{ind},2); 
                ind = find(~cellfun('isempty',B.coeffs),1);   
                nrB = size(B.coeffs{ind},1); 
                ncB = size(B.coeffs{ind},2); 
                if ncA ~= nrB
                    error('dimensions of matrices must agree')
                end   
                KNc = get_KNg(A.vertices,A.degree+B.degree); 
                KNb = get_KNg(B.vertices,B.degree);
                for i=1:get_JNg(A.vertices,A.degree+B.degree)
                    Cmon = zeros(nrA,ncB);
                    j = KNc(i,:);
                    for m=1:get_JNg(B.vertices,B.degree)
                        k = KNb(m,:);
                        if compk(j,k) % if j>=k elementwise    
                            Cmon = Cmon + eval(['A.coeffs{' gMS(j-k) '}*B.coeffs{' gMS(k) '};']);
                        end
                    end
                    eval(['coeffs{' gMS(j) '} = Cmon;']);
                end
                C = A.cl(coeffs);                
            end
        end
        
        function Ad = increase_degree(A,d,vertices)
            if d == 0 % if no increase in degree
                Ad = A;
                return
            end
%             if isa(A,'double') % if A is a constant matrix
%                 if nargin < 3
%                     error('number of vertices unknown');
%                 end
%                 KNd = get_KNg(vertices,d);
%                 for i=1:get_JNg(vertices,d)
%                     j  = KNd(i,:);    
%                     eval(['coeffs{' gMS(j) '} = (factorial(d)/Pi(j))*A;']);
%                 end
%             else % increase degree of A(a) to g+d
                ind  = find(~cellfun('isempty',A.coeffs),1);
                nr   = size(A.coeffs{ind},1); 
                nc   = size(A.coeffs{ind},2);
                KNgd = get_KNg(A.vertices,A.degree+d); 
                KNd  = get_KNg(A.vertices,d);
                for i=1:get_JNg(A.vertices,A.degree+d)
                    Amon = zeros(nr,nc);
                    j = KNgd(i,:);    
                    for m=1:get_JNg(A.vertices,d)
                        k = KNd(m,:);        
                        if compk(j,k) % if j>=k elementwise
                            Amon = Amon + (factorial(d)/Pi(k))*eval(['A.coeffs{' gMS(j-k) '};']);
                        end
                    end
                    eval(['coeffs{' gMS(j) '} = Amon;'])
                end
%             end
            Ad = A.cl(coeffs);
        end
    end
end