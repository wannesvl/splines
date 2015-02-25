function p = parameter(n)
    % Creates n linear polynomial functions in n variables
    if n == 1
        p = Polynomial({0, 1}');
        return
    end
    data = zeros(n, 2^n);
    data(:, 1 + 2.^(0:n-1)) = eye(n)
    coeffs = Coefficients(data, [1, 2^n], [n, 1]);
    I = reshape(1:2^n, 2 * ones(1, n));
    coeffs = coeffs(I)
    p = Polynomial(coeffs);
