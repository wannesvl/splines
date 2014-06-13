function monomial_string = get_monomial_string(k) 
% get monomial string from vector
k = k+1;
monomial_string = num2str(k(1));
if length(k)>1
    for i = 2:length(k)
        monomial_string = [ monomial_string ',' num2str(k(i)) ];
    end
end
end