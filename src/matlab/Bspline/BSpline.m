classdef BSpline
    properties (Access=protected)
        cl
    end
    properties
        coeffs
        basis
    end
    methods
        function s = BSpline(basis, coeffs)
            % Constructor for BSpline
            % validateattributes(coeffs(:), {'double'}, {'size', [length(basis), 1]}, 'BSpline', 'coeffs');
            s.basis = basis;
            if isa(coeffs, 'BSplineCoeffs')
                s.coeffs = coeffs;
            else
                s.coeffs = BSplineCoeffs(coeffs);
            end
            % Validate input
            % if length(s.basis) == 1 && length(s.coeffs) ~= length(s.basis)
            % if any(size(s.coeffs.coeffs) ~= arrayfun(@length, s.basis))
            %     error('Coefficient array is not of correct size')
            % end
            % Idea: make coeffs a seperate class
            % to handle vector and matrix valued splines
            s.cl = str2func(class(s));
        end

        function c = getcoeffs(self)
            c = self.coeffs.coeffs;
            if self.coeffs.isscalar || self.coeffs.isvector
                c = vertcat(c{:});
            end
        end

        function s = f(self, x)
            % Evaluate B-spline at x
            c = self.basis.f(x) * self.coeffs;
            s = c.coeffs;
            if self.coeffs.isscalar
                s = vertcat(s{:});
            end
        end

        function s = plus(self, other)
            if isa(self, class(other))
                basis = self.basis + other.basis;
                coeffs = basis.transform(self.basis) * self.coeffs + ...
                         basis.transform(other.basis) * other.coeffs;
                s = self.cl(basis, coeffs);
            else
                try
                    basis = self.basis;
                    coeffs = self.coeffs + other;
                    s = self.cl(basis, coeffs);
                catch err
                    s = other + self;
                end
            end
        end

        function s = uminus(self)
            s = self.cl(self.basis, -self.coeffs);
        end

        function s = minus(self, other)
            s = self + (- other);
        end

        function s = mtimes(self, other)
            if isa(self, class(other))
                basis = self.basis * other.basis;
                grev = basis.greville;
                b_self = self.basis.f(grev);
                b_other = other.basis.f(grev);
                [i, j] = self.basis.pairs(other.basis);
                basis_product = b_self(:, i) .* ...
                                b_other(:, j);
                T = basis.f(grev) \ basis_product;
                T(abs(T) < 1e-10) = 0;
                coeffs_product = self.coeffs(i) .* other.coeffs(j);
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

        function s = vertcat(varargin)
            % Concatenate two splines

            % Determine support of concatenation
            k_min = nan;
            k_max = nan;
            for i=1:length(varargin)
                if isa(varargin{i}, 'BSpline')
                    k_min = min(k_min, varargin{i}.basis.knots(1));
                    k_max = max(k_max, varargin{i}.basis.knots(end));
                end
            end
            % Convert constants to Bspline bases
            for i=1:length(varargin)
                if isa(varargin{i}, 'double')
                    varargin{i} = BSpline(BSplineBasis([k_min, k_max], 0), {varargin{i}}); %#ok<CCAT1>
                end
            end
            % Now concatenate the bases
            b = varargin{1}.basis;
            for i=2:length(varargin)
                b = varargin{i}.basis + b;
            end
            c = cell(size(varargin));
            for i=1:length(varargin)
                c{i} = b.transform(varargin{i}.basis) * varargin{i}.coeffs;
            end
            c = vertcat(c{:});
            s = varargin{1}.cl(b, c);
        end

        function s = horzcat(varargin)
            % Concatenate two splines

            % Determine support of concatenation
            k_min = nan;
            k_max = nan;
            for i=1:length(varargin)
                if isa(varargin{i}, 'BSpline')
                    k_min = min(k_min, varargin{i}.basis.knots(1));
                    k_max = max(k_max, varargin{i}.basis.knots(end));
                end
            end
            % Convert constants to Bspline bases
            for i=1:length(varargin)
                if isa(varargin{i}, 'double')
                    varargin{i} = BSpline(BSplineBasis([k_min, k_max], 0), {varargin{i}});
                end
            end
            % Now concatenate the bases
            b = varargin{1}.basis;
            for i=2:length(varargin)
                b = varargin{i}.basis + b;
            end
            c = cell(size(varargin));
            for i=1:length(varargin)
                c{i} = b.transform(varargin{i}.basis) * varargin{i}.coeffs;
            end
            c = horzcat(c{:});
            s = varargin{1}.cl(b, c);
        end

        function s = transpose(self)
            s = self.cl(self.basis, self.coeffs');
        end

        function s = ctranspose(self)
            s = self.cl(self.basis, self.coeffs.');
        end

        function d = derivative(self, o)
            if nargin == 1
                o = 1;
            end
            B, P = self.basis.derivative(o);
            d = self.cl(B, P * self.coeffs);
        end

        function s = insert_knots(self, knots)
            basis = self.basis.insert_knots(knots);
            coeffs = basis.transform(self.basis) * self.coeffs;
            s = self.cl(basis, coeffs);
        end

        function s = increase_degree(self, degree)
            basis = self.basis.increase_degree(degree);
            coeffs = basis.transform(self.basis) * self.coeffs;
            s = self.cl(basis, coeffs);
        end

        function i = integral(self)
            k = self.basis.knots;
            d = self.basis.degree;
            i = (k(d + 2:end) - k(1:end - d - 1))'  / (d + 1) * self.coeffs;
            i = i.coeffs{1};
        end

        function t = trace(self)
            % Sum of the diagonal coeffs
            for i=1:length(self.basis)
                c(i) = trace(self.coeffs.coeffs{i});
            end
            t = self.cl(self.basis, c);
        end
    end
end