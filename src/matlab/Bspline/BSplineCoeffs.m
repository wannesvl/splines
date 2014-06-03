classdef BSplineCoeffs
    properties (Access=protected)
        cl
        shape
    end
    properties
        coeffs
    end
    methods
        function c = BSplineCoeffs(coeffs)
            % Constructor for BSplineCoeffs
            %
            % The coefficients are either vector for scalar valued
            % splines, or cells of vectors and matrices for vector
            % or matrix valued splines
            %
            % This is not the preferred way. Better to make class for
            % scalar, vector and matrix valued coefficients???
            if isa(coeffs, 'cell')
                if isvector(coeffs)
                    coeffs = coeffs(:);
                end
                % Check if all elements are of equal size
                sizes = cellfun(@size, coeffs, 'UniformOutput', false);
                if length(sizes) | isequal(sizes{:})
                    c.coeffs = coeffs;
                else
                    error('Coefficients should all be of equal size');
                end
            elseif isvector(coeffs)  % Scalar valued
                coeffs = coeffs(:);
                c.coeffs = mat2cell(coeffs, ones(size(coeffs, 1), 1));
            else
                error('Coeffs not correctly formatted');
            end
            c.shape = size(c.coeffs{1});
            c.cl = str2func(class(c));
        end

        function l = length(self)
            l = length(self.coeffs);
        end

        function s = size(self)
            s = size(self.coeffs);
        end

        function b = isscalar(self)
            % Are we dealing with scalar coefficients?
            b = isscalar(self.coeffs{1});
        end

        function b = isvector(self)
            % Are we dealing with vector coefficients?
            b = isvector(self.coeffs{1});
        end

        function c = plus(self, other)
            % Returns the sum of coefficients
            %
            % Returns: 
            %    BSplineCoeffs: The sum of self.coeffs and other.coeffs.
            %       If either of the terms is a matrix the function sums each of
            %       the coefficient with the matrix
            if isa(self, class(other))
                c = self.cl(cellfun(@plus, self.coeffs, other.coeffs, 'UniformOutput', false));
            else
                try
                    c = self.cl(cellfun(@(v) v + other, self.coeffs, 'UniformOutput', false));
                catch err
                    c = other + self;
                end
            end
        end

        function c = sum(self)
            % Only for univariate splines
            c = ones(length(self), 1)' * self;
            c = c.coeffs{1};
        end

        function c = uminus(self)
            c = self.cl(cellfun(@uminus, self.coeffs, 'UniformOutput', false));
        end

        function c = times(self, other)
            % pointwise product
            if isa(self, class(other))
                c = cellfun(@mtimes, self.coeffs, other.coeffs, 'UniformOutput', false);
            else
                try
                    c = cellfun(@(v) v * other, self.coeffs, 'UniformOutput', false);
                catch err
                    c = cellfun(@(v) self * v, other.coeffs, 'UniformOutput', false);
                end
            end
            c = BSplineCoeffs(c);
        end

        function c = mtimes(a, self)
            % Implement 'generalized' inner product
            if isa(a, 'double')
                if isvector(self.coeffs)
                    a = kron(a, eye(self.shape(1)));
                    % s = self.size;
                    % a = repmat(a, prod(s(2:end)));
                    c = a * vertcat(self.coeffs{:});
                    % And now convert back to cell
                    c = mat2cell(c, ... 
                        self.shape(1) * ones(size(c, 1) / self.shape(1), 1), size(c, 2));
                    c = self.cl(c);
                else
                    a = kron(a, eye(self.shape(1)))
                    s = self.size;
                    % a = repmat(a, prod(s(2:end)));
                    c = a * reshape(vertcat(self.coeffs{:}), [], s(2) * self.shape(2));
                    % And now convert back to cell
                    c = mat2cell(c, ... 
                        self.shape(1) * ones(size(c, 1) / self.shape(1), 1), size(c, 2));
                    c = self.cl(c);
                end
            else  % Recursive implementation
                c = 0;
            end
        end

        function c = vertcat(varargin)
            % Concatenate matrices vertically
            d = cellfun(@(v) v.coeffs, varargin, 'UniformOutput', false);
            c = BSplineCoeffs(cellfun(@vertcat, d{:}, 'UniformOutput', false));
        end

        function c = horzcat(varargin)
            % Concatenate matrices horizontaly
            d = cellfun(@(v) v.coeffs, varargin, 'UniformOutput', false);
            c = BSplineCoeffs(cellfun(@horzcat, d{:}, 'UniformOutput', false));
        end

        function c = transpose(self)
            % Transpose each of the coefficients
            %
            % What about multivariate spline coefficients??
            c = self.cl(cellfun(@transpose, self.coeffs, 'UniformOutput', false)');
        end

        function c = ctranspose(self)
            % Conjugate transpose each of the coefficients
            c = self.cl(cellfun(@ctranspose, self.coeffs, 'UniformOutput', false)');
        end

        function varargout = subsref(self, s)
            % Note to this function: c.coeffs{:} does not work!!
            % Instead you should define dummy d = c.coeffs and d{:}
            %
            % For this reason I consider changing the implementation to object arrays
            if strcmp(s(1).type, '.')
                if any(strcmp(s(1).subs, properties(self))) || ...
                   any(strcmp(s(1).subs, methods(self)))
                    [varargout{1:nargout}] = builtin('subsref', self, s);
                else
                    error(['''%s'' is not a public property or method ' ...
                        'of the BSplineCoeffs class.'],s(1).subs);
                end
            else
                varargout{1} = self.cl(subsref(self.coeffs, s));
            end
        end
    end
end