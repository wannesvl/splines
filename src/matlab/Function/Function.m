classdef Function
    properties (Access=protected)
        cl
    end
    properties
        coeffs
        basis
    end
    methods
        function s = Function(basis, coeffs)
            % A generic function class
            % 
            % Function(b, c) creates a function object defined by bases b and
            % coefficients c.
            %
            % Args:
            %    basis (cell of basis objects): The function bases to use     
            %    coeffs (array, cell or Coefficients): The function
            %      coefficients. Array input is allowed only for scalar coefficients
            %
            % Returns:
            %    An instance of the function class
            if isa(basis, 'UnivariateBasis') && isscalar(basis)
                basis = {basis};
            end
            s.basis = basis;
            lengths = cellfun(@length, s.basis);
            if isa(coeffs, 'Coefficients')
                s.coeffs = coeffs;
            else
                if all(size(coeffs) == lengths) && ~isa(coeffs, 'cell') % Scalar coefficients
                    % sp = arrayfun(@(i) ones(size(coeffs, i), 1), 1:length(basis), 'UniformOutput', false);
                    % coeffs = mat2cell(coeffs, sp{:});
                    coeffs = num2cell(coeffs);
                end
                s.coeffs = Coefficients(coeffs);
            end
            % Validate input
            if lengths ~= size(s.coeffs)
                error('Function coefficient of different size than bases')
            end
            s.cl = str2func(class(s));
        end

        function s = f(self, x)
            % Evaluate a Function at x
            warning('OFF', 'MATLAB:mat2cell:TrailingUnityVectorArgRemoved')
            if self.dims == 1 && ~isa(x, 'cell')
                s = self.basis{1}.f(x) * self.coeffs;
            else
                s = cellfun(@(b, x) b.f(x), self.basis, x, 'UniformOutput', false) * self.coeffs;
            end
            if self.coeffs.isscalar  % If scalar coefficients, convert to regular matrix
                s = s.coeffs2tensor;
            else  % Return cell array
                s = s.coeffs;
            end
        end

        function s = f_partial(self, x, idx)
            % Partial evaluation of coordinate idx of f
            %
            % Returns a new Function object
            b = arrayfun(@(i) self.basis{i}.f(x{i}), idx, 'UniformOutput', false);
            
            s = 0;
        end


        function d = dims(self)
            d = length(self.basis);
        end

        function s = plus(self, other)
            if isa(self, class(other))
                basis = cellfun(@plus, self.basis, other.basis, 'UniformOutput', false);
                Tself = cellfun(@(b1, b2) b1.transform(b2), basis, self.basis, 'UniformOutput', false);
                Tother = cellfun(@(b1, b2) b1.transform(b2), basis, other.basis, 'UniformOutput', false);
                coeffs = Tself * self.coeffs + Tother * other.coeffs;
                s = self.cl(basis, coeffs);
            else  % This implementation depends on the type of basis we are dealing with and cannot be implemented
                error('plus with other than Function objects is not defined. Consider defining plus in a subclass.')
            end
        end

        function s = uminus(self)
            s = self.cl(self.basis, -self.coeffs);
        end

        function s = minus(self, other)
            s = self + (- other);
        end

        function s = mtimes(self, other)
            % This is a general implementation but could be performed more efficiently in subclasses
            % It is recommended to overload mtimes in the subclass
            function T = transform(A, B)
                T = A \ B;
                T(abs(T) < 1e-10) = 0;
            end
            if isa(self, class(other)) 
                % Basis of product
                basis = cellfun(@mtimes, self.basis, other.basis, 'UniformOutput', false);
                % Take kronecker product of coefficients
                [i_other, i_self] = arrayfun(@(i) find(ones(size(other.coeffs, i), size(self.coeffs, i))), 1:self.dims, 'UniformOutput', false);  % Give all indices of products
                coeffs_product = self.coeffs(i_self{:}) * other.coeffs(i_other{:});
                % Determine transformation matrices
                x = cellfun(@(b) b.x_, basis, 'UniformOutput', false);
                b = cellfun(@(b, x) b.f(x), basis, x, 'UniformOutput', false);
                b_self = cellfun(@(b, x) b.f(x), self.basis, x, 'UniformOutput', false);
                b_other = cellfun(@(b, x) b.f(x), other.basis, x, 'UniformOutput', false);
                basis_product = cellfun(@(b1, b2, is, io) b1(:, is) .* b2(:, io), b_self, b_other, i_self, i_other, 'UniformOutput', false);
                T = cellfun(@(b, bi) transform(b, bi), b, basis_product, 'UniformOutput', false);
                s = self.cl(basis, T * coeffs_product);
            else   % Assume multiplication with array
                try
                    basis = self.basis;
                    coeffs = self.coeffs .* other;
                    s = self.cl(basis, coeffs);
                catch err
                    basis = other.basis;
                    coeffs = self .* other.coeffs;
                    s = other.cl(basis, coeffs);
                end
            end
        end

        % function s = times(self, other)
        %     % Elementwise multiplication for vector or matrix valued splines
        %     function T = transform(A, B)
        %         T = A \ B;
        %         T(abs(T) < 1e-10) = 0;
        %     end
        %     if any(self.coeffs.shape ~= other.coeffs.shape)
        %         error('Coefficient shape are not compatible')
        %     end
        %     basis = cellfun(@mtimes, self.basis, other.basis, 'UniformOutput', false);
        %     % Take kronecker product of coefficients
        %     [i_other, i_self] = arrayfun(@(i) find(ones(size(other.coeffs, i), size(self.coeffs, i))), 1:self.dims, 'UniformOutput', false);  % Give all indices of products
        %     coeffs_product = self.coeffs(i_self{:}) .* other.coeffs(i_other{:});
        %     % Determine transformation matrices
        %     x = cellfun(@(b) b.x_, basis, 'UniformOutput', false);
        %     b = cellfun(@(b, x) b.f(x), basis, x, 'UniformOutput', false);
        %     b_self = cellfun(@(b, x) b.f(x), self.basis, x, 'UniformOutput', false);
        %     b_other = cellfun(@(b, x) b.f(x), other.basis, x, 'UniformOutput', false);
        %     basis_product = cellfun(@(b1, b2, is, io) b1(:, is) .* b2(:, io), b_self, b_other, i_self, i_other, 'UniformOutput', false);
        %     T = cellfun(@(b, bi) transform(b, bi), b, basis_product, 'UniformOutput', false);
        %     s = self.cl(basis, T * coeffs_product);
        % end

        function s = horzcat(varargin)
            l = cellfun(@(b) length(b.basis), varargin, 'ErrorHandler', @(varargin) NaN);
            l = l(~isnan(l));
            if ~all(l == l(1))
                error('Cannot concatenate: the basis have different size')
            end

            % Find a common basis for all terms
            b = cell(1, l(1));
            for i=1:length(varargin)
                if isa(varargin{i}, mfilename)
                    b = cellfun(@plus, b, varargin{i}.basis, 'UniformOutput', false);
                    func = str2func(class(varargin{i}));
                end
            end
            size_b = cellfun(@length, b);
            if isscalar(size_b)  % Correction for univariate splines
                size_b = [size_b, 1];
            end
            
            % Compute coefficients in common basis
            c = cell(size(varargin));
            for i=1:length(varargin)
                if isa(varargin{i}, mfilename)
                    T = cellfun(@(b, bi) b.transform(bi), b, varargin{i}.basis, 'UniformOutput', false);
                    c{i} = T * varargin{i}.coeffs;
                else  % Constant function: Simply repeat matrices along dimensions of b
                    c{i} = Coefficients(repmat({varargin{i}}, size_b));
                end
            end

            % Finally concatenate all coefficients
            c = horzcat(c{:});
            s = func(b, c);
        end

        function s = vertcat(varargin)
            l = cellfun(@(b) length(b.basis), varargin, 'ErrorHandler', @(varargin) NaN);
            l = l(~isnan(l));
            if ~all(l == l(1))
                error('Cannot concatenate: the basis have different size')
            end

            % Find a common basis for all terms
            b = cell(1, l(1));
            for i=1:length(varargin)
                if isa(varargin{i}, mfilename)
                    b = cellfun(@plus, b, varargin{i}.basis, 'UniformOutput', false);
                    func = str2func(class(varargin{i}));
                end
            end
            size_b = cellfun(@length, b);
            if isscalar(size_b)  % Correction for univariate splines
                size_b = [size_b, 1];
            end
            
            % Compute coefficients in common basis
            c = cell(size(varargin));
            for i=1:length(varargin)
                if isa(varargin{i}, mfilename)
                    T = cellfun(@(b, bi) b.transform(bi), b, varargin{i}.basis, 'UniformOutput', false);
                    c{i} = T * varargin{i}.coeffs;
                else  % Constant function: Simply repeat matrices along dimensions of b
                    c{i} = Coefficients(repmat({varargin{i}}, size_b));
                end
            end

            % Finally concatenate all coefficients
            c = vertcat(c{:});
            s = func(b, c);
        end

        function s = transpose(self)
            s = self.cl(self.basis, self.coeffs.');
        end

        function s = ctranspose(self)
            s = self.cl(self.basis, self.coeffs');
        end

        function varargout = subsref(self, s)
            if strcmp(s(1).type, '.')
                if any(strcmp(s(1).subs, properties(self))) || ...
                   any(strcmp(s(1).subs, methods(self)))
                    [varargout{1:nargout}] = builtin('subsref', self, s);
                else
                    error(['''%s'' is not a public property or method ' ...
                        'of the Coefficients class.'], s(1).subs);
                end
            elseif strcmp(s(1).type, '()')
                basis = self.basis;
                coeffs = cellfun(@(c) builtin('subsref', c, s(1)), self.coeffs.coeffs, 'UniformOutput', false);
                if isscalar(s)
                    varargout{1} = self.cl(basis, Coefficients(coeffs));
                else
                    [varargout{1:nargout}] = builtin('subsref', self.cl(basis, Coefficients(coeffs)), s(2:end));
                end
            else
                error('Invalid use of {}')
            end
        end

        function last = end(self, idx, n)
            % Overload end
            last = length(self.coeffs.shape(idx));
        end

        function b = ge(self, other)
            % Overload comparisons to allow more intuitive constraints
            %
            % Only works for YALMIP at the moment
            b = [];
            if isa(self, 'double')
                temp = self;
                self = other;
                other = temp;
            end
            c = self.coeffs.coeffs;
            if ~isvector(c{1})
                for i=1:numel(c)
                    b = [b, 0.5 * (c{i} + c{i}') >= other];
                end
            else
                for i=1:numel(c)
                    b = [b, c{i} >= other];
                end
            end
        end

        function b = le(self, other)
            % Overload comparisons to allow more intuitive constraints
            %
            % Only works for YALMIP at the moment
            b = [];
            if isa(self, 'double')
                temp = self;
                self = other;
                other = temp;
            end
            c = self.coeffs.coeffs;
            if ~isvector(c{1})
                for i=1:numel(c)
                    b = [b, 0.5 * (c{i} + c{i}') <= other];
                end
            else
                for i=1:numel(c)
                    b = [b, c{i} <= other];
                end
            end
        end

        function b = eq(self, other)
            % Overload comparisons to allow more intuitive constraints
            %
            % Only works for YALMIP at the moment
            b = [];
            if isa(self, 'double')
                temp = self;
                self = other;
                other = temp;
            end
            c = self.coeffs.coeffs;
            for i=1:numel(c)
                b = [b, c{i} == other];
            end
        end

        function s = double(self)
            % Overload double for use with YALMIP variables
            s = self.cl(self.basis, double(self.coeffs));
        end

        function s = increase_degree(self, degree)
            % Increase each basis by degree(i)
            %
            % Args: 
            %    degree (vector): Degree increase for each basis
            %
            % Returns:            
            %    BSpline: Bspline with updated coefficients
            b = arrayfun(@(i) self.basis{i}.increase_degree(degree(i)), 1:self.dims, 'UniformOutput', false);
            T = cellfun(@(b1, b2) b1.transform(b2), b, self.basis, 'UniformOutput', false);
            s = self.cl(b, T * self.coeffs);
        end

        function s = trace(self)
            % Return the trace of a matrix valued spline
            c = cellfun(@trace, self.coeffs.coeffs, 'UniformOutput', false);
            s = self.cl(self.basis, c);
        end
    end

    methods (Static)
        function s = sdpvar(basis, dim, p)
            if isscalar(basis) && ~isa(basis, 'cell')
                basis = {basis};
            end
            dim = num2cell(dim);
            cl = class(basis{1});
            cl = str2func(cl(1:end-5));  % Determine type of function (Polynomial, BSpline, ...)
            lengths = cellfun(@length, basis);
            if isscalar(lengths)
                lengths = [lengths, 1];
            end
            if nargin == 2
                coeffs = arrayfun(@(i) sdpvar(dim{:}), zeros(lengths), 'UniformOutput', false);
            elseif nargin == 3
                coeffs = arrayfun(@(i) sdpvar(dim{:}, p), zeros(lengths), 'UniformOutput', false);
            end
            s = cl(basis, coeffs);
        end
    end
end
