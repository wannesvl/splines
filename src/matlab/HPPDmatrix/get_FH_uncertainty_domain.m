function [F,H] = get_FH_uncertainty_domain(N,b,type)

if strcmp(type,'6v') % use polytopic domain with 6 vertices
    switch b
        case 0 
            V = [ eye(N); zeros(N,N) ];
        case 1
            V = zeros(2*N,N*N);
            for i = 1:1:N
                V(i,(i-1)*N+1) = 1;
                ind = 1;
                for j = 1:1:N
                    if j ~= i
                        V([i N+i N+j],(i-1)*N+ind+1) = [1; -1; 1];
                        ind = ind+1;
                    end
                end
            end
        otherwise
            V = zeros(2*N,N*N+N*(N-1));
            for i = 1:1:N
                V(i,(i-1)*N+1) = 1;
                k = 1;
                for j = 1:1:N
                    if j ~= i
                        V([i N+i N+j],(i-1)*N+k+1) = [1; -b; b];
                        V([i j N+i N+j],N^2+(i-1)*(N-1)+k) = [b; 1-b; -b; b];
                        k = k+1;
                    end
                end
            end
    end
    F = V(1:N,:);
    H = V(N+1:end,:);
elseif strcmp(type,'3v') % divide polytopic domain in subdomains with 3 vertices
    % MANUALLY IMPLEMENTED NOW FOR SPECIFIC SUBDIVISION...
    if N > 2
        error('only works for N = 2 vertices')
    end
    switch b
        case 0
            [F,H] = deal(cell(1,1));
            F{1} = eye(2);
            H{1} = zeros(2,2);
        case 1
            [F,H] = deal(cell(1,2));
            F{1} = [0 0 1; 1 1 0];
            F{2} = [0 1 1; 1 0 0];
            H{1} = [0 1 0; 0 -1 0];
            H{2} = [0 0 -1; 0 0 1];
        otherwise
            [F,H] = deal(cell(1,4));
            F{1} = [0 0 1-b; 1 1 b];
            F{2} = [0 1-b 1; 1 b 0];
            F{3} = [0 1 b; 1 0 1-b];
            F{4} = [1 1 b; 0 0 1-b];
            H{1} = [0 b b; 0 -b -b];
            H{2} = [0 b 0; 0 -b 0];
            H{3} = [0 0 -b; 0 0 b];
            H{4} = [0 -b -b; 0 b b];
    end   
else 
    error('incorrect (subdivision) of polytopic domain')
end