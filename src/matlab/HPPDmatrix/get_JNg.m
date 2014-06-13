function JNg = get_JNg(N,g) % cardinality of KNg
    JNg = factorial(N+g-1)/(factorial(g)*factorial(N-1));
end