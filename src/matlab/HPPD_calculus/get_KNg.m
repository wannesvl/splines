function [KNg] = get_KNg(N,g)

%generates all the solutions of: k_1 + k_2 + ... + k_N = g
KNg = [];
if N == 1
    KNg = g;
else
    for i = 0:g
        temp = get_KNg(N-1,g-i);
        KNg = [ KNg
            i*ones(size(temp,1),1) temp ];
    end
end
