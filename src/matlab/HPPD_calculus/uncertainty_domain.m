function [F,H] = uncertainty_domain(N,b)

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


