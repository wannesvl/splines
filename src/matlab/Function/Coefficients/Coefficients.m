classdef Coefficients
    properties (Access=protected)
        cl
    end
    % properties (Access={?Coefficients,?Function})
    %     shape
    % end
    properties
        coeffs
        shape
    end
    methods
        function c = Coefficients(coeffs)
            % Constructor for Coefficients
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

        function s = size(self, i)
            if nargin == 1
                s = size(self.coeffs);
            else
                s = size(self.coeffs, i);
            end
        end

        function b = isscalar(self)
            % Are we dealing with scalar coefficients?
            b = isscalar(self.coeffs{1});
        end

        function b = isvector(self)
            % Are we dealing with vector coefficients?
            b = isvector(self.coeffs{1});
        end

        function m = coeffs2tensor(self)
            % Convert the coefficients to tensor
            %
            % Shortened copy of mat2cell
            c = self.coeffs;
            elements = numel(c);

            if elements == 1
                m = c{1};
                return
            end

            csize = size(c);
            % Construct the matrix by concatenating each dimension of the cell array into
            %   a temporary cell array, CT
            % The exterior loop iterates one time less than the number of dimensions,
            %   and the final dimension (dimension 1) concatenation occurs after the loops

            % Loop through the cell array dimensions in reverse order to perform the
            %   sequential concatenations
            for cdim=(length(csize)-1):-1:1
                % Pre-calculated outside the next loop for efficiency
                ct = cell([csize(1:cdim) 1]);
                cts = size(ct);
                ctsl = length(cts);
                mref = {};

                % Concatenate the dimension, (CDIM+1), at each element in the temporary cell
                %   array, CT
                for mind=1:prod(cts)
                    [mref{1:ctsl}] = ind2sub(cts,mind);
                    % Treat a size [N 1] array as size [N], since this is how the indices
                    %   are found to calculate CT
                    if ctsl==2 && cts(2)==1
                        mref = {mref{1}};
                    end
                    % Perform the concatenation along the (CDIM+1) dimension
                    ct{mref{:}} = cat(cdim+1,c{mref{:},:});
                end
                % Replace M with the new temporarily concatenated cell array, CT
                c = ct;
            end
            % Finally, concatenate the final rows of cells into a matrix
            m = cat(1,c{:});
        end

        function c = plus(self, other)
            % Returns the sum of coefficients
            %
            % Returns: 
            %    Coefficients: The sum of self.coeffs and other.coeffs.
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

        function c = minus(self, other)
            c = self + (-other);
        end

        function c = times(self, other)
            % pointwise product
            if isa(self, class(other))
                c = cellfun(@times, self.coeffs, other.coeffs, 'UniformOutput', false);
            else
                if isa(self, 'Coefficients')
                    c = cellfun(@(v) v * other, self.coeffs, 'UniformOutput', false);
                else
                    c = cellfun(@(v) self * v, other.coeffs, 'UniformOutput', false);
                end
            end
            c = Coefficients(c);
        end

        function c = mtimes(A, self)
            % Matrix multiplication with same type or transformation matrices
            if isa(A, class(self))  % Both inputs are coefficients
                c = self.cl(cellfun(@mtimes, A.coeffs, self.coeffs, 'UniformOutput', false));
                return
            end
            coeffs = self.coeffs2tensor;
            if isa(A, 'double')  % Univariate transformation matrix
                if isvector(self.coeffs)
                    coeffs = kron(A, eye(self.shape(1))) * coeffs;
                else
                    error('Univariate splines cannot be multiplied with cell array of matrices')
                end
            elseif isa(A, 'cell')  % Multivariate transformation matrices
                if isvector(self.coeffs)
                    A = cellfun(@(a) kron(a, eye(self.shape(1))), A, 'UniformOutput', false);
                else  % Not yet correct for dims >= 3
                    A = cellfun(@(a, i) kron(a, eye(self.shape(i))), A, num2cell(1:ndims(self.coeffs)), 'UniformOutput', false, 'ErrorHandler', @(varargin) varargin{2});
                end
                if length(A) == 1
                    coeffs = A{1} * coeffs;
                else
                    coeffs = tmprod(coeffs, A, 1:ndims(self.coeffs));
                end
            end

            % The tricky part is converting it to a cell again
            % D = cell(1, ndims(self.coeffs));
            D = {};
            for i=1:ndims(self.coeffs)
                if i > 2
                    D{i} = ones(size(coeffs, i), 1);
                else
                    D{i} = ones(size(coeffs, i) / self.shape(i), 1);
                end
            end
            D{1} = self.shape(1) * D{1};
            D{2} = self.shape(2) * D{2};
            c = self.cl(mat2cell(coeffs, D{:}));
        end

        function c = multiply(self, A, dims)
            % Multiply coeffs with cell array A along dimensions dims
            %
            % 
            coeffs = self.coeffs2tensor;
            size1 = cellfun(@(a) size(a, 1), A)
            A = cellfun(@(a) kron(a, eye(self.shape(1))), A, 'UniformOutput', false);
            coeffs = squeeze(tmprod(coeffs, A, dims))
            s = size(self);
            c = 0;
            return
            if ndims(self.coeffs) - length(dims) > 0
                reshape(coeffs, s(1:end - dims))
            end
            c = 0;
            % D = {};
            % for i=1:ndims(self.coeffs)
            %     D{i} = ones(size(coeffs, i))
        end

        function c = vertcat(varargin)
            % Concatenate matrices vertically
            d = cellfun(@(v) v.coeffs, varargin, ...
                        'UniformOutput', false, ...
                        'ErrorHandler', @(s, i) num2cell(i));  % Concatenation with arrays
            c = Coefficients(cellfun(@vertcat, d{:}, 'UniformOutput', false));
        end

        function c = horzcat(varargin)
            % Concatenate matrices horizontaly
            d = cellfun(@(v) v.coeffs, varargin, ...
                        'UniformOutput', false, ...
                        'ErrorHandler', @(s, i) num2cell(i));
            c = Coefficients(cellfun(@horzcat, d{:}, 'UniformOutput', false));
        end

        function c = transpose(self)
            % Transpose each of the coefficients
            %
            % What about multivariate spline coefficients??
            c = self.cl(cellfun(@transpose, self.coeffs, 'UniformOutput', false));
        end

        function c = ctranspose(self)
            % Conjugate transpose each of the coefficients
            c = self.cl(cellfun(@ctranspose, self.coeffs, 'UniformOutput', false));
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
                        'of the Coefficients class.'],s(1).subs);
                end
            else
                varargout{1} = self.cl(subsref(self.coeffs, s));
            end
        end

        function c = double(self)
            c = self.cl(cellfun(@double, self.coeffs, 'UniformOutput', false));
        end
    end
end