classdef Coefficients
    properties
        data
        siz
        shape
    end

    properties (Access=protected)
        cl
    end

    methods
        function blktens = Coefficients(varargin)
            if nargin == 3
                data = varargin{1};
                blktens.siz = varargin{2};
                blktens.shape = varargin{3};
                if ~strcmp(data, 'empty')
                    blktens.data = reshape(data, blktens.siz(1) * blktens.shape(1), []);
                end
            elseif nargin == 2
                % infer size from shape
                data = varargin{1};
                sdata = size(data);
                shape = varargin{2};
                blktens.data = reshape(data, size(data, 1), []);
                blktens.shape = shape;
                blktens.siz = [sdata(1) / shape(1), sdata(2) / shape(2), sdata(3:end)];
            end
            blktens.cl = str2func(mfilename);
        end

        function tens = astensor(self)
            tens = reshape(self.data, self.totalsize);
        end

        function siz = size(self, i)
            if nargin == 1
                siz = self.siz;
            else
                siz = self.siz(i);
            end
        end


        function siz = totalsize(self)
            siz = [self.siz(1) * self.shape(1), ...
                   self.siz(2) * self.shape(2), ...
                   self.siz(3:end)];
        end

        function bool = isscalar(self)
            bool = all(self.shape == [1, 1]);
        end

        function bool = isvector(self)
            bool = self.shape(2) == 1;
        end

        function blktens = plus(self, other)
            if isa(self, class(other))
                blktens = self.cl(self.data + other.data, self.siz, self.shape);
            elseif isa(self, mfilename)
                other = repmat(other, [self.siz(1), prod(self.size(2:end))]);
                blktens = self.cl(self.data + other, self.siz, self.shape);
            else
                blktens = other + self;
            end
        end

        function blktens = times(self, other)
            % Multiply each block of self with the corresponding block from other
            %
            % [A B; C D] .* [E F; G H] = [AE BF; CG DH]

            % Creat all zero tensor and populate using linear indexing
            if isa(self, class(other))
                blktens = self.cl(zeros(size(self.data)), self.siz, self.shape);
                for i = 1:prod(self.siz)
                    % For some reason we cannot use overloaded method so we
                    % have to call subsref directly
                    S = struct('type', {'()', '.'}, 'subs', {{i}, 'data'});
                    blktens = blktens.subsasgn(S(1), self.subsref(S) .* other.subsref(S));
                end
            elseif isa(self, mfilename)
                blktens = self.cl(zeros(size(self.data)), self.siz, self.shape);
                for i = 1:prod(self.siz)
                    S = struct('type', {'()', '.'}, 'subs', {{i}, 'data'});
                    blktens = blktens.subsasgn(S(1), self.subsref(S) * other);
                end
            else
                blktens = other.cl(zeros(size(other.data)), other.siz, other.shape);
                for i = 1:prod(other.siz)
                    S = struct('type', {'()', '.'}, 'subs', {{i}, 'data'});
                    blktens = blktens.subsasgn(S(1), self * other.subsref(S));
                end
            end
        end

        function blktens = mtimes(T, other)
            % Full mode tensor matrix product
            if isa(T, mfilename)
                self = T;
                % This is not correct!
                blktens = self.cl(zeros(size(self.data)), self.siz, self.shape);
                [m, n] = size(self.data);
                blktens = self.cl(sdpvar(m,n), self.siz, self.shape);
                for i = 1:prod(self.siz)
                    S = struct('type', {'()', '.'}, 'subs', {{i}, 'data'});
                    blktens = blktens.subsasgn(S(1), self.subsref(S) .* other.subsref(S));
                end
                return
            end
            if ~iscell(T)
                T = {T};
            end
            blktens = other.tmprod(T, 1:length(T));
        end

        function blktens = tmprod(self, U, mode)
            % Inspired by tmprod in Tensorlab.

            % kron matrices with I to account for block shape
            for i=1:min(length(U), 2)
                U{i} = kron(U{i}, speye(self.shape(i)));
            end

            % Heuristically sort modes
            size_tens = self.totalsize;
            [~, idx] = sort(size_tens(mode) ./ cellfun('size', U, 1));
            mode = mode(idx);
            U = U(idx);

            % Prepermute tensor
            n = length(mode);
            N = length(self.siz);
            modec = setxor(mode, 1:N);
            perm = [mode modec];
            size_tens = size_tens(perm);
            S = self.data;
            I = reshape(1:prod(size_tens), size_tens);
            if any(mode ~= 1:n)
                I = reshape(permute(I, perm), size_tens(1), []);
                S = S(I);
            end

            % Cycle through the products
            for i = 1:n
                if ~isscalar(U{i})
                    size_tens(1) = size(U{i}, 1);
                end
                I = reshape(1:prod(size_tens), size_tens);
                I = permute(I, [2:N, 1]);
                S = U{i} * S;
                if i < n
                    S = S(reshape(I, size(S, 1), []));
                    size_tens = size_tens([2:N 1]);
                end
            end

            % Inverse permute
            iperm(perm([n:N 1:n-1])) = 1:N;
            siz = [size_tens(1:2) ./ self.shape, size_tens(3:end)];
            I = reshape(1:prod(size_tens), size_tens(iperm));
            I = permute(I, iperm);
            S = S(reshape(I, size(S, 1), []));
            blktens = self.cl(S, siz(iperm), self.shape);
        end

        function blktens = vertcat(varargin)
            % Find out which input is of blktens class
            for i = 1:length(varargin)
                if isa(varargin{i}, mfilename)
                    t = varargin{i};
                    break
                end
            end
            S = struct('type', {'()', '.'}, 'subs', {{':'}, 'data'});
            data = [];
            len = 0;
            for i = 1:length(varargin)
                if isa(varargin{i}, mfilename)
                    data = [data; varargin{i}.subsref(S)];
                    len = len + varargin{i}.shape(1);
                else
                    data = [data; repmat(varargin{i}, 1, prod(t.siz))]
                    len = len + size(varargin{i}, 1);
                end
            end
            blktens = t.cl(zeros(size(data)), t.siz, [len, t.shape(2)]);
            S = struct('type', {'()'}, 'subs', {{':'}});
            blktens = blktens.subsasgn(S, data);
        end

        function blktens = horzcat(varargin)
            % Find out which input is of blktens class
            for i = 1:length(varargin)
                if isa(varargin{i}, mfilename)
                    t = varargin{i};
                    break
                end
            end
            S = struct('type', {'()', '.'}, 'subs', {{':'}, 'data'});
            data = [];
            len = 0;
            for i = 1:length(varargin)
                if isa(varargin{i}, mfilename)
                    data = [data varargin{i}.subsref(S)'];
                    len = len + varargin{i}.shape(2);
                else
                    data = [data repmat(varargin{i}, prod(t.siz), 1)]
                    len = len + size(varargin{i}, 2);
                end
            end
            data'
            blktens = t.cl(zeros(size(data)), t.siz, [t.shape(1), len]);
            S = struct('type', {'()'}, 'subs', {{':'}});
            blktens = blktens.subsasgn(S, data');
        end

        function blktens = ctranspose(self)
            % Transpose each of the blocks separately
            blktens = self.cl(zeros(size(self.data)), self.siz, [self.shape(2), self.shape(1)]);
            for i = 1:prod(self.siz)
                S = struct('type', {'()', '.'}, 'subs', {{i}, 'data'});
                blktens = blktens.subsasgn(S(1), self.subsref(S)');
            end
        end

        function n = end(self)
            n = prod(self.siz)
        end

        function varargout = subsref(self, S)
            switch S(1).type
                case '()'
                    if length(S(1).subs) == 1  % Linear indexing
                        if strcmp(S(1).subs{1}, ':')
                            S(1).subs{1} = 1:prod(self.siz);
                        end
                        [i, j] = ind2sub([self.siz(1), prod(self.size(2:end))], S(1).subs{1});
                        I = repmat((1:self.shape(1))', size(i, 1), size(i, 2) * self.shape(2));
                        J = repmat(1:self.shape(2), size(j, 1) * self.shape(1), size(j, 2));
                        i = kron(i, ones(self.shape));
                        j = kron(j, ones(self.shape));
                        idx = sub2ind(size(self.data), (i - 1) * self.shape(1) + I, (j - 1) * self.shape(2) + J);
                        y = self.cl(self.data(idx), self.shape);
                    else % Convert indices to linear indices?
                        idx = sub2ind(self.siz, S(1).subs{:});
                        S(1).subs = {idx};
                        [varargout{1:nargout}] = subsref(self, S);
                        return
                    end
                    if length(S) > 1
                        [varargout{1:nargout}] = subsref(y, S(2:end));
                    else
                        varargout{1} = y;
                    end
                case '.'
                    if any(strcmp(S(1).subs, properties(self))) || ...
                       any(strcmp(S(1).subs, methods(self)))
                        [varargout{1:nargout}] = builtin('subsref', self, S);
                    end
            end
        end

        function y = subsasgn(self, S, varargin)
            switch S(1).type
                case '()'
                    if length(S(1).subs) == 1  % Linear indexing
                        if strcmp(S(1).subs{1}, ':')
                            S(1).subs{1} = 1:prod(self.siz);
                        end
                        [i, j] = ind2sub([self.siz(1), prod(self.size(2:end))], S(1).subs{1});
                        I = repmat((1:self.shape(1))', size(i, 1), size(i, 2) * self.shape(2));
                        J = repmat(1:self.shape(2), size(j, 1) * self.shape(1), size(j, 2));
                        i = kron(i, ones(self.shape));
                        j = kron(j, ones(self.shape));
                        idx = sub2ind([self.shape .* self.siz(1:2), self.siz(3:end)], (i - 1) * self.shape(1) + I, (j - 1) * self.shape(2) + J);
                        y = self;
                        y.data(idx) = varargin{1};
                    else
                        error('Not yet implemented')
                    end
                case '.'
                    if any(strcmp(S(1).subs, properties(self))) || ...
                       any(strcmp(S(1).subs, methods(self)))
                        y = builtin('subsasgn', self, S, varargin);
                    end
            end
        end
    end
end
