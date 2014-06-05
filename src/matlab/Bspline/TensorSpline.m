classdef TensorSpline
    properties (Access=protected)
        cl
    end
    properties
        coeffs
        basis
    end
    methods
        function s = TensorSpline(basis, coeffs)
            % Constructor for TensorSpline
            s.basis = basis;
            lengths = cellfun(@length, s.basis);
            if isa(coeffs, 'BSplineCoeffs')
                s.coeffs = coeffs;
            else
                if size(coeffs) == lengths % Scalar coefficients
                    coeffs = mat2cell(coeffs, ones(size(coeffs, 1), 1), ones(size(coeffs, 2), 1));
                end
                s.coeffs = BSplineCoeffs(coeffs);
            end
            % Validate input
            if lengths ~= size(s.coeffs)
                error('B-spline coefficient of different size than basis')
            end
            s.cl = str2func(class(s));
        end

        function s = f(self, x)
            s = cellfun(@(b, x) b.f(x), self.basis, x, 'UniformOutput', false) * self.coeffs;
            if self.coeffs.isscalar
                s = s.coeffs2tensor;
            end
        end

        function d = dims(self)
            d = length(self.basis);
        end

        function s = plus(self, other)
            if isa(other, 'TensorSpline')
                basis = cellfun(@plus, self.basis, other.basis, 'UniformOutput', false);
                Tself = cellfun(@(b1, b2) b1.transform(b2), basis, self.basis, 'UniformOutput', false);
                Tother = cellfun(@(b1, b2) b1.transform(b2), basis, other.basis, 'UniformOutput', false);
                coeffs = Tself * self.coeffs + Tother * other.coeffs;
            else
                basis = self.basis;
                coeffs = self.coeffs + other;
            end
            s = self.cl(basis, coeffs);
        end

        function s = mtimes(self, other)
            function T = transform(A, B)
                T = A \ B;
                T(abs(T) < 1e-10) = 0;
            end
            if isa(self, class(other))
                basis = cellfun(@mtimes, self.basis, other.basis, 'UniformOutput', false);
                [p_self, p_other] = cellfun(@(b1, b2) b1.pairs(b2), self.basis, other.basis, 'UniformOutput', false);
                grev = cellfun(@(b) b.greville, basis, 'UniformOutput', false);
                b = cellfun(@(b, g) b.f(g), basis, grev, 'UniformOutput', false);
                b_self = cellfun(@(b, g) b.f(g), self.basis, grev, 'UniformOutput', false);
                b_other = cellfun(@(b, g) b.f(g), other.basis, grev, 'UniformOutput', false);
                basis_product = cellfun(@(b1, b2, ps, po) b1(:, ps) .* b2(:, po), b_self, b_other, p_self, p_other, 'UniformOutput', false);
                coeffs_product = self.coeffs(p_self{:}) .* other.coeffs(p_other{:});
                T = cellfun(@(b, bi) transform(b, bi), b, basis_product, 'UniformOutput', false);
                s = self.cl(basis, T * coeffs_product);
            else
                try
                    basis = self.basis;
                    coeffs = self.coeffs .* other;
                    s = self.cl(basis, coeffs);
                catch err
                    % s = other * self
                    basis = other.basis;
                    coeffs = self .* other.coeffs;
                    s = other.cl(basis, coeffs);
                end
            end
        end

        function s = uminus(self)
            s = self.cl(self.basis, -self.coeffs);
        end

        function s = minus(self, other)
            s = self + (- other);
        end

        function s = horzcat(varargin)
            l = cellfun(@(b) length(b.basis), varargin, 'ErrorHandler', @(varargin) NaN);
            l = l(~isnan(l));
            if ~all(l == l(1))
                error('Cannot concatenate: the basis have different size')
            end

            % Find a common basis for all terms
            b = cell(1, l);
            for i=1:length(varargin)
                if isa(varargin{i}, mfilename)
                    b = cellfun(@plus, b, varargin{i}.basis, 'UniformOutput', false);
                end
            end
            size_b = cellfun(@length, b);
            
            % Compute coefficients in common basis
            c = cell(size(varargin));
            for i=1:length(varargin)
                if isa(varargin{i}, mfilename)
                    T = cellfun(@(b, bi) b.transform(bi), b, varargin{i}.basis, 'UniformOutput', false);
                    c{i} = T * varargin{i}.coeffs;
                else  % Constant function: Simply repeat matrices along dimensions of b
                    c{i} = BSplineCoeffs(repmat({varargin{i}}, size_b));
                end
            end

            % Finally concatenate all coefficients
            c = horzcat(c{:});
            cl = str2func(mfilename);
            s = cl(b, c);
        end

        function s = vertcat(varargin)
            l = cellfun(@(b) length(b.basis), varargin, 'ErrorHandler', @(varargin) NaN);
            l = l(~isnan(l));
            if ~all(l == l(1))
                error('Cannot concatenate: the basis have different size')
            end

            % Find a common basis for all terms
            b = cell(1, l);
            for i=1:length(varargin)
                if isa(varargin{i}, mfilename)
                    b = cellfun(@plus, b, varargin{i}.basis, 'UniformOutput', false);
                end
            end
            size_b = cellfun(@length, b);
            
            % Compute coefficients in common basis
            c = cell(size(varargin));
            for i=1:length(varargin)
                if isa(varargin{i}, mfilename)
                    T = cellfun(@(b, bi) b.transform(bi), b, varargin{i}.basis, 'UniformOutput', false);
                    c{i} = T * varargin{i}.coeffs;
                else  % Constant function: Simply repeat matrices along dimensions of b
                    c{i} = BSplineCoeffs(repmat({varargin{i}}, size_b));
                end
            end

            % Finally concatenate all coefficients
            c = vertcat(c{:});
            cl = str2func(mfilename);
            s = cl(b, c);
        end

        function s = transpose(self)
            s = self.cl(self.basis, self.coeffs');
        end

        function s = ctranspose(self)
            s = self.cl(self.basis, self.coeffs.');
        end

        function i = integral(self)
            i = cellfun(@(b) (b.knots(b.degree + 2:end) - b.knots(1:end - b.degree - 1))' ...
                        / (b.degree + 1), self.basis, 'UniformOutput', false) * self.coeffs;
            i = i.coeffs{1};
        end
    end
end
