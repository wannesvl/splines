classdef BSplineCoeffs
    properties (Access=protected)
        cl
        size
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
            %
            % TODO: validate and "tidy" incorrect input
            %
            % Currently only works 
            if isa(coeffs, 'cell')
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
            c.size = size(c.coeffs{1});
            c.cl = str2func(class(c));
        end

        function l = length(self)
            l = length(self.coeffs);
        end

        function b = isscalar(self)
            % Are we dealing with scalar coefficients?
            b = isscalar(self.coeffs{1});
        end

        function b = isvector(self)
            % Are we dealing with scalar coefficients?
            b = isvector(self.coeffs{1});
        end

        function c = plus(self, other)
            % Sum of two coefficients
            if isa(self, class(other))
                % c = vertcat(self.coeffs{:}) + vertcat(other.coeffs{:});
                % c = mat2cell(c, self.size(1) * ones(length(c) / self.size(1), 1), size(c, 2));
                % c = self.cl(reshape(c, size(self.coeffs)));
                c = self.cl(cellfun(@plus, self.coeffs, other.coeffs, 'UniformOutput', false));
            else
                try
                    % c = vertcat(self.coeffs{:}) + other;
                    % c = mat2cell(c, self.size(1) * ones(length(c) / self.size(1), 1), size(c, 2));
                    % c = self.cl(reshape(c, size(self.coeffs)));
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
                    c = kron(a, eye(self.size(1))) * vertcat(self.coeffs{:});
                    % And now convert back to cell
                    c = mat2cell(c, ... 
                        self.size(1) * ones(size(c, 1) / self.size(1), 1), size(c, 2));
                    c = self.cl(c);
                end
            else  % Recursive implementation

            end
        end

        function c = vertcat(varargin)
            % Concatenate matrices
            d = cellfun(@(v) v.coeffs, varargin, 'UniformOutput', false);
            c = BSplineCoeffs(cellfun(@vertcat, d{:}, 'uni', false));
        end

        function c = horzcat(varargin)
            % Concatenate matrices
            d = cellfun(@(v) v.coeffs, varargin, 'UniformOutput', false);
            c = BSplineCoeffs(cellfun(@horzcat, d{:}, 'uni', false));
        end

        function c = transpose(self)
            % Transpose each of the coefficients
            c = self.cl(cellfun(@transpose, self.coeffs, 'UniformOutput', false));
        end

        function c = ctranspose(self)
            % Transpose each of the coefficients
            c = self.cl(cellfun(@ctranspose, self.coeffs, 'UniformOutput', false));
        end

        function varargout = subsref(self, s)
            % Note to this function: c.coeffs{:} does not work!!
            % Instead you should define dummy d = c.coeffs and d{:}
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