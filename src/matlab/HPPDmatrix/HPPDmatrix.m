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
            % Also works for constant matrices:
            % -> coeffs is a double
            % -> degree = 0
            % -> no vertices
            %
            % For HPPD matrices: 
            % -> parameter in unit simplex
            % -> coeffs is a cell array where the element
            %    coeffs{k_1+1,...,k_N+1} corresponds to monomial 
            %    k = (k_1,...,k_N), where k_1 + ... + k_N = vertices.
            %    All other cell elements are assumed to be empty.
            % -> vertices and degree are derived from coeffs.
            %
            % Ex.: A(a1,a2) = a1 A_10 + a2 A_01
            % is stored in a 2x2 cell as A{2,1} = A_10, A{1,2} = A_01,
            % and represents a HPPD matrix of degree 1 with 2 vertices
            if isa(coeffs,'double')   % constant matrix
                A.degree   = 0;
                A.coeffs   = coeffs;
                A.cl       = str2func(class(A));
            elseif isa(coeffs,'cell') % HPPD matrix
                A.vertices = length(size(coeffs));
                A.degree   = size(coeffs,1) - 1;
                A.coeffs   = coeffs;
                A.cl       = str2func(class(A));
            else
                error('inconsistent class for coeffs')
            end
        end
        
        function b = ge(k1,k2) % check if k1 >= k2 elementwise
            if length(k1)==length(k2)
                b=0;
                if sum(k1>=k2)==length(k1),
                    b=1;
                end
            else
                error('Cannot compare monomial coefficients of different length!')
            end
        end
             
        function s = plus(A,B)
            if A.vertices ~= B.vertices
                error('self and other should have same number of vertices!')
            end
            if A.degree ~= B.degree % homogenize s.t. A and B have same degree
                if A.degree > B.degree && B.degree > 0
                   B = increase_degree(B,A.degree-B.degree); 
                elseif A.degree > B.degree && B.degree == 0
                   B = increase_degree(B,A.degree-B.degree,A.vertices);
                elseif B.degree > A.degree && A.degree > 0
                   A = increase_degree(A,B.degree-A.degree);
                else
                   A = increase_degree(A,B.degree-A.degree,B.vertices); 
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
            if isa(A.coeffs,'double') % constant matrix
                coeffs = -A.coeffs;
            else                      % HPPD matrix
                coeffs = cellfun(@uminus, A.coeffs, 'UniformOutput', false);
            end
            C = A.cl(coeffs);
        end      
        
        function C = minus(A,B)
            C = A + (-B);
        end
        
        function B = ctranspose(A)
            if isa(A.coeffs,'double') % constant matrix
                coeffs = (A.coeffs)';
            else                      % HPPD matrix
                coeffs = cellfun(@ctranspose, A.coeffs, 'UniformOutput', false);
            end
            B = A.cl(coeffs);
        end      
        
        function C = mtimes(A,B)       
            % The case that A and/or B is a constant matrix is taken into account  
            if isa(A.coeffs,'double') && isa(B.coeffs,'double')
                try 
                    coeffs = A.coeffs*B.coeffs;
                    C = A.cl(coeffs);
                catch err 
                    coeffs = A.coeffs.*B.coeffs;
                end
            elseif isa(A.coeffs,'double') && isa(B.coeffs,'cell')
                ind = find(~cellfun('isempty',B.coeffs),1);   
                KNg = get_KNg(B.vertices,B.degree);
                if size(A.coeffs,2) ~= size(B.coeffs{ind},1)
                    try              
                        for j = 1:get_JNg(B.vertices,B.degree)
                            jj = get_monomial_string(KNg(j,:));
                            eval(['coeffs{' jj '} = A.coeffs.*B.coeffs{' jj '};'])
                        end
                    catch err
                        error('dimensions of matrices must agree')
                    end
                else             
                    for j = 1:get_JNg(B.vertices,B.degree)
                        jj = get_monomial_string(KNg(j,:));
                        eval(['coeffs{' jj '} = A.coeffs*B.coeffs{' jj '};'])
                    end
                end 
                C = B.cl(coeffs); 
            elseif isa(A.coeffs,'cell') && isa(B.coeffs,'double')
                ind = find(~cellfun('isempty',A.coeffs),1);  
                KNg = get_KNg(A.vertices,A.degree);                
                if size(A.coeffs{ind},2) ~= size(B.coeffs,1)
                    try              
                        for j = 1:get_JNg(B.vertices,B.degree)
                            jj = get_monomial_string(KNg(j,:));
                            eval(['coeffs{' jj '} = A.coeffs{' jj '}.*B.coeffs;'])
                        end
                    catch err
                        error('dimensions of matrices must agree')
                    end
                else
                    for j = 1:get_JNg(A.vertices,A.degree)
                        jj = get_monomial_string(KNg(j,:));
                        eval(['coeffs{' jj '} = A.coeffs{' jj '}*B.coeffs;'])
                    end
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
                        if j >= k % if j>=k elementwise    
                            Cmon = Cmon + eval(['A.coeffs{' get_monomial_string(j-k) '}*B.coeffs{' get_monomial_string(k) '};']);
                        end
                    end
                    eval(['coeffs{' get_monomial_string(j) '} = Cmon;']);
                end
                C = A.cl(coeffs);                
            end
        end
        
        function Ad = increase_degree(A,d,vertices)
            if d == 0 % if no increase in degree
                Ad = A;
                return
            end
            if isa(A.coeffs,'double') % if A is a constant matrix
                if nargin < 3
                    error('number of vertices unknown');
                end
                KNd = get_KNg(vertices,d);
                for i=1:get_JNg(vertices,d)
                    j  = KNd(i,:);    
                    eval(['coeffs{' get_monomial_string(j) '} = (factorial(d)/pi_coeff(j))*A.coeffs;']);
                end
            else % increase degree of A(a) to vertices+d
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
                        if j >= k % if j>=k elementwise
                            Amon = Amon + (factorial(d)/pi_coeff(k))*eval(['A.coeffs{' get_monomial_string(j-k) '};']);
                        end
                    end
                    eval(['coeffs{' get_monomial_string(j) '} = Amon;'])
                end
            end
            Ad = A.cl(coeffs);
        end
        
        function s = horzcat(varargin)
            % horizontal concatenation of HPPD matrices
            
            % determine maximum polynomial degree of all terms
            for i = 1:length(varargin)
               if isa(varargin{i}, 'HPPDmatrix')
                   degrees(i) = varargin{i}.degree;
               else 
                   error('inconsistent class for input argument')
               end
            end
            [max_degree, index] = max(degrees);
            vertices = varargin{index}.vertices;
            
            % homogenize all terms 
            for i = 1:length(varargin)
                varargin{i} = varargin{i}.increase_degree(max_degree-degrees(i),vertices);
            end
            
            % concatenation
            dimension_coeffs = (varargin{1}.degree + 1)*ones(1,vertices);
            coeffs = cell(1,prod(dimension_coeffs));
            if max_degree == 0
                term = [];
                for i = 1:length(varargin)
                    term = [term, varargin{i}.coeffs];
                end
                coeffs = term;
                s = varargin{1}.cl(coeffs);
            else
                index = find(~cellfun('isempty',varargin{1}.coeffs));
                for j = 1:length(index)
                    term = [];
                    for i = 1:length(varargin)
                        term = [term, varargin{i}.coeffs{index(j)}];
                    end
                    coeffs{index(j)} = term;
                end
                coeffs = reshape(coeffs,dimension_coeffs);
                s = varargin{1}.cl(coeffs);
            end
        end
        
        function s = vertcat(varargin)
            % vertical concatenation of HPPD matrices
            
            % determine maximum polynomial degree of all terms
            for i = 1:length(varargin)
               if isa(varargin{i}, 'HPPDmatrix')
                   degrees(i) = varargin{i}.degree;
               else 
                   error('inconsistent class for input argument')
               end
            end
            [max_degree, index] = max(degrees);
            vertices = varargin{index}.vertices;
            
            % homogenize all terms 
            for i = 1:length(varargin)
                varargin{i} = varargin{i}.increase_degree(max_degree-degrees(i),vertices);
            end
                  
            % concatenation
            dimension_coeffs = (varargin{1}.degree + 1)*ones(1,vertices);
            coeffs = cell(1,prod(dimension_coeffs));
            if max_degree == 0
                term = [];
                for i = 1:length(varargin)
                    term = [term; varargin{i}.coeffs];
                end
                coeffs = term;
                s = varargin{1}.cl(coeffs);
            else
                index = find(~cellfun('isempty',varargin{1}.coeffs));
                for j = 1:length(index)
                    term = [];
                    for i = 1:length(varargin)
                        term = [term; varargin{i}.coeffs{index(j)}];
                    end
                    coeffs{index(j)} = term;
                end
                coeffs = reshape(coeffs,dimension_coeffs);
                s = varargin{1}.cl(coeffs);
            end               
        end
        
        function c = getcoeffs(A)
            index = find(~cellfun('isempty',A.coeffs));
            c = cell(1,length(index));
            for j = 1:length(index)
                c{j} = A.coeffs{index(j)};
            end
        end
        
        function s = f(A,x) % evaluate HPPD matrix at x (in unit simplex)
            KNg = get_KNg(A.vertices,A.degree);
            s = 0;
            for j = 1:length(A.getcoeffs);
               jj = KNg(j,:);
               s = s + eval(['prod(x.^jj)*A.coeffs{' get_monomial_string(jj) '};']);
            end
        end      
        
        function s = trace(A)
            KNg = get_KNg(A.vertices,A.degree);
            for j = 1:length(A.getcoeffs)
                jj = get_monomial_string(KNg(j,:));
                eval(['coeffs{' jj '} = trace(A.coeffs{' jj '});']);
            end
            s = A.cl(coeffs);
        end
        
        function s = integral(A)
            s = 0;
            for j = 1:get_JNg(A.vertices,A.degree)
                 s = s + A.getcoeffs{j}/nchoosek(get_JNg(A.vertices,A.degree)-1,j-1);
            end
            s = s / (A.degree + 1);
        end
            
        % TO DO: 
        % - evaluation of HPPD matrix at different points simultaneously
        % - transformation to gamma-domain for given bound b (0 < b < 1) 
        %   on rate of parameter variation
        %   1) 6 vertices
        %   2) 4x 3 vertices
    end
end