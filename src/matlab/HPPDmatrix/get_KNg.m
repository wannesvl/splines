function KNg = get_KNg(N,g)
%generate all the solutions of: k1 + k2 + ... + kN = g
KNg = [];
if N == 1
    KNg = g;
else
    for i = 0:g
        temp = get_KNg(N-1,g-i);
        KNg = [ KNg; i*ones(size(temp,1),1), temp ];
    end
end
end