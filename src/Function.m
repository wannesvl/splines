
% splines - A framework for manipulating and optimizing (B-)splines
% Copyright (C) 2015 Wannes Van Loock
%
% This file is part of splines.
%
% splines is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% splines is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU Lesser General Public License for more details.
%
% You should have received a copy of the GNU Lesser General Public License
% along with splines.  If not, see <http://www.gnu.org/licenses/>.

classdef (InferiorClasses = {?casadi.MX,?casadi.SX}) Function
    properties (Constant)
        BLKDIV = 100;
    end
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
                    % coeffs = num2cell(coeffs);
                    siz = size(coeffs);
                    shape = [1,1];
                    data = coeffs;
                elseif isa(coeffs, 'cell')
                    siz = size(coeffs);
                    shape = size(coeffs{1});
                    data = cell2mat(coeffs);
                end
                s.coeffs = Coefficients(data, siz, shape);
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
            if self.dims == 1 && ~isa(x, 'cell') && ~isa(x, 'Function')
                s = self.basis{1}.f(x) * self.coeffs;
            elseif isa(self, 'Polynomial') && (isa(x, 'Function') || isa(x{1}, 'Function'))
                % evaluate scalar-valued polynomial at a function
                s = 0;
                basiseval = cellfun(@(b, x) b.f(x), self.basis, x, 'UniformOutput', false);
                sbs = cell(length(self.coeffs.siz), 1);
                for i=1:prod(self.coeffs.siz)
                    if self.coeffs(i).data ~= 0
                        [sbs{:}] = ind2sub(self.coeffs.siz, i);
                        p = basiseval{1}{sbs{1}};
                        for j=2:length(sbs)
                            p = p * basiseval{j}{sbs{j}};
                        end
                        s = s + p * self.coeffs(i).data;
                    end
                end
                return
            else
                s = cellfun(@(b, x) b.f(x), self.basis, x, 'UniformOutput', false) * self.coeffs;
            end
            if self.coeffs.isscalar %&& ~isa(s,'BSpline') % If scalar coefficients, convert to regular matrix
                s = s.astensor;
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
                    % varargin{i}
                    c{i} = Coefficients(repmat(varargin{i}, size_b), size_b, size(varargin{i}));
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
                elseif isempty(varargin{i})
                    break
                else  % Constant function: Simply repeat matrices along dimensions of b
                    c{i} = Coefficients(repmat(varargin{i}, size_b), size_b, size(varargin{i}));
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

        function s = diag(self, v)
            % Create a matrix valued spline from a vector valued one
            if nargin == 1
                v = 0;
            end
            s = self.cl(self.basis, diag(self.coeffs, v));
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
                prodsiz = prod(self.coeffs.siz);
                prodshape = prod(self.coeffs.shape);
                c = self.coeffs(1:prodsiz).data;
                if length(s(1).subs) == 2
                    idx = sub2ind(self.coeffs.shape, s(1).subs{:});
                else
                    idx = s(1).subs{1};
                end
                % size(repmat(idx, prodsiz, 1)), size(kron((0:prodshape:prodshape*prodsiz)', ones(size(idx))))
                I = repmat(idx, prodsiz, 1) + kron((0:prodshape:prodshape*prodsiz-1)', ones(size(idx)));
                coeffs = Coefficients(c(I), [prodsiz, 1], size(idx));
                I = reshape(1:prodsiz, self.coeffs.siz);
                coeffs = coeffs(I);
                %bla
                % Needs vectorization!!
                % coeffs = builtin('subsref', self.coeffs, s(1));
                % c1 = builtin('subsref', self.coeffs(1).data, s(1));
                % coeffs = Coefficients(repmat(c1, 1, prod(self.coeffs.siz)), self.coeffs.siz, size(c1));
                % for i=1:prod(self.coeffs.siz)
                %     c = self.coeffs(i).data;
                %     coeffs(i) = builtin('subsref', c, s(1));
                % end
                % coeffs = cellfun(@(c) builtin('subsref', c, s(1)), self.coeffs.coeffs, 'UniformOutput', false);
                if isscalar(s)
                    varargout{1} = self.cl(basis, coeffs);
                else
                    [varargout{1:nargout}] = builtin('subsref', self.cl(basis, coeffs), s(2:end));
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
            if isa(self, 'numeric')
                b = le(other, self);
                return
            end
            c = self.coeffs;
            if isa(c.data, 'sdpvar') || isa(c.data, 'numeric')
                if ~isvector(c(1))
                    nel = prod(c.siz);
                    b = [];
                    for i = 1:self.BLKDIV:nel
                        m = min(i+self.BLKDIV-1, nel);
                        d = c(i:m).spblkdiag();
                        b = [b, 0.5 * (d + d') >= other];
                    end
                    b = unblkdiag(b);
                    % for i=1:prod(c.siz)
                    %     b = [b, 0.5 * (c(i).data + c(i).data') >= other];
                    % end
                else
                    b = [c(:).data >= repmat(other, prod(c.totalsize) / numel(other), 1)];
                    % for i=1:prod(c.siz)
                    %     b = [b, c(i).data >= other];
                    % end
                end
            else
                b = Constraint(c.data(:), other, inf);
            end
        end

        function b = le(self, other)
            % Overload comparisons to allow more intuitive constraints
            %
            % Only works for YALMIP at the moment
            b = [];
            if isa(self, 'double')
                b = ge(other, self);
                return
            end
            c = self.coeffs;
            if isa(c.data, 'sdpvar') || isa(c.data, 'numeric')
                if ~isvector(c(1))
                    nel = prod(c.siz);
                    b = [];
                    for i = 1:self.BLKDIV:nel
                        m = min(i+self.BLKDIV, nel);
                        d = c(i:m).spblkdiag();
                        b = [b, 0.5 * (d + d') <= other];
                    end
                    b = unblkdiag(b);
                    % for i=1:prod(c.siz)
                    %     b = [b, 0.5 * (c(i).data + c(i).data') <= other];
                    % end
                else
                    b = [c(:).data <= repmat(other, prod(c.totalsize) / numel(other), 1)];
                    % for i=1:prod(c.siz)
                    %     b = [b, c(i).data <= other];
                    % end
                end
            else
                b = Constraint(c.data(:), -inf, other);
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
            c = self.coeffs;
            if isa(c.data, 'sdpvar') || isa(c.data, 'numeric')
                for i=1:prod(c.siz)
                    b = [b, c(i).data == other];
                end
            else
                b = Constraint(c.data(:), other, other);
            end
        end

        function s = value(self)
            % Overload double for use with YALMIP variables
            s = self.cl(self.basis, value(self.coeffs));
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
            % c = cellfun(@trace, self.coeffs.coeffs, 'UniformOutput', false);
            s = self.cl(self.basis, trace(self.coeffs));
        end
    end

    methods (Static)
        function s = sdpvar(basis, shape, p)
            if isscalar(basis) && ~isa(basis, 'cell')
                basis = {basis};
            end
            % dim = num2cell(dim);
            cl = class(basis{1});
            cl = str2func(cl(1:end-5));  % Determine type of function (Polynomial, BSpline, ...)
            lengths = cellfun(@length, basis);
            if isscalar(lengths)
                lengths = [lengths, 1];
            end
            if nargin == 2
                coeffs = sdpvar(shape(1) * lengths(1), shape(2) * prod(lengths(2:end)));
                coeffs = Coefficients(coeffs, lengths, shape);
            elseif nargin == 3  % This is incorrect! should be applied on subblocks!
                coeffs = sdpvar(shape(1) * ones(1, prod(lengths)), shape(2) * ones(1, prod(lengths)), p);
                coeffs = Coefficients(horzcat(coeffs{:}), [1, prod(lengths)], shape);
                I = reshape(1:prod(lengths), lengths);
                coeffs = coeffs(I);
            end
            s = cl(basis, coeffs);
        end

        function s = MX(basis, shape, name)
            if isscalar(basis) && ~isa(basis, 'cell')
                basis = {basis};
            end
            % dim = num2cell(dim);
            cl = class(basis{1});
            cl = str2func(cl(1:end-5));  % Determine type of function (Polynomial, BSpline, ...)
            lengths = cellfun(@length, basis);
            if isscalar(lengths)
                lengths = [lengths, 1];
            end
            coeffs = casadi.MX.sym(name, shape(1) * lengths(1), shape(2) * prod(lengths(2:end)));
            coeffs = Coefficients(coeffs, lengths, shape);
            s = cl(basis, coeffs);
        end
    end
end
