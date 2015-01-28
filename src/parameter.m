function p = parameter(varargin)
    % Create the polynomial f(x) = x along each dimension
    if length(varargin) == 1
        dim = varargin{1};
        if length(dim) == 1
            m = dim; n = 1;
        else
            m = dim(1); n = dim(2);
        end
    elseif length(varargin) == 2
        m = varargin{1}; n = varargin{2};
    else
        error('Invalid parameter dimension')
    end

    switched = false;
    if m < n
        temp = m; m = n; n = temp;
        switched = true;
    end
    if m * n > 1
        p_coeffs = cell(2 * ones(1, m * n));
        [p_coeffs{:}] = deal(zeros(m, n));
        P = num2cell(eye(m * n) + 1);
        count = 1;
        for i=1:m
            for j=1:n
                temp = zeros(m, n);
                temp(i, j) = 1;
                p_coeffs{P{:, count}} = temp;
                count = count + 1;
            end
        end
    else
        p_coeffs = {0, 1};
    end
    p = Polynomial(Coefficients(cell2mat(p_coeffs), size(p_coeffs), size(p_coeffs{1})));
    if switched
        p = p';
    end
end
