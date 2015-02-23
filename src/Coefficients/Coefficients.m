% TODO: replace subsref/subsasgn loops by vectorized versions!

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
                blktens.data = reshape(data, blktens.siz(1) * blktens.shape(1), []);
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

       % function blktens = set.siz(self, value)
       %     blktens = self.cl(self.data, value, self.shape);
       % end

        function tens = astensor(self)
            tens = reshape(self.data, self.totalsize);
        end

        function blktens = ascell(self)
            % Return cell representation of self
            blktens = mat2cell(self.data, self.shape(1) * ones(1, self.siz(1)), self.shape(2) * ones(1, prod(self.siz(2:end))));
            blktens = reshape(blktens, self.siz);
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
            bool = self.shape(2) == 1 || self.shape(1) == 1;
        end

        % function blktens = reshape(self, siz)
        %     % Reshape the tensor according to size parameter
        %     I = reshape(1:prod(self.siz), siz)
        %     blktens = blktens.subsref(struct('type', {'()'}, 'subs', {{I}}));
        % end

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

        function blktens = uminus(self)
            blktens = self.cl(-self.data, self.siz, self.shape);
        end

        function blktens = minus(self, other)
            blktens = self + (-other);
        end

        function blktens = times(self, other)
            % Multiply each block of self with the corresponding block from other
            if isa(self, class(other))
                blktens = self.cl(self.data .* other.data, self.siz, self.shape);
            elseif isa(self, mfilename)
                if isscalar(other)  % Speedup for scalar multiplication
                    blktens = other.cl(self.data * other, self.siz, self.shape);
                    return
                end
                if isscalar(self)  % Correct size
                    data = kron(self.spblkdiag(), eye(size(other, 1))) * repmat(other, prod(self.siz), 1);
                    blktens = self.cl(data, [prod(self.siz), 1], size(other));
                else
                    data = self.spblkdiag() * repmat(other, prod(self.siz), 1);
                    blktens = self.cl(data, [prod(self.siz), 1], [self.shape(1), size(other, 2)]);
                end
                I = reshape(1:prod(self.siz), [self.siz(1), prod(self.siz(2:end))]);
                blktens = blktens.subsref(struct('type', {'()'}, 'subs', {{I}}));
                blktens.siz = self.siz;
            else
                if isscalar(self)
                    blktens = other.cl(self * other.data, other.siz, other.shape);
                    return
                end
                if isscalar(other)
                    data = repmat(self, 1, prod(other.siz)) * kron(other.spblkdiag(), eye(size(self, 1)));
                    blktens = other.cl(data, [1, prod(other.siz)], size(self));
                else
                    data = repmat(self, 1, prod(other.siz)) * other.spblkdiag();
                    blktens = other.cl(data, [1, prod(other.siz)], [size(self, 1), other.shape(2)]);
                end
                I = reshape(1:prod(other.siz), [other.siz(1), prod(other.siz(2:end))]);
                blktens = blktens.subsref(struct('type', {'()'}, 'subs', {{I}}));
                blktens.siz = other.siz;
            end
        end

        function blktens = mtimes(T, other)
            % Full mode tensor matrix product
            % Code requires thorough checks!!!!
            if isa(T, mfilename)
                % [A B; C D] .* [E F; G H] = [AE BF; CG DH]
                self = T;
                % This is not correct!
                % blktens = self.cl(zeros(size(self.data)), self.siz, self.shape);
                % [m, n] = size(self.data);
                % blktens = self.cl(sdpvar(m,n), self.siz, self.shape);
                % blktens = other;
                % for i = 1:prod(self.siz)
                %     S = struct('type', {'()', '.'}, 'subs', {{i}, 'data'});
                %     blktens = blktens.subsasgn(S(1), self.subsref(S) * other.subsref(S));
                % end
                % return
                S = struct('type', {'()', '.'}, 'subs', {{':'}, 'data'});
                % other = other';  % Is this necessary?
                if isnumeric(self.data)
                    data = self.spblkdiag() * other.subsref(S);
                    blktens = self.cl(data, [prod(self.siz), 1], [self.shape(1), other.shape(2)]);
                else
                    data = self.subsref(S)' * other.spblkdiag()';
                    blktens = self.cl(data, [1, prod(self.siz)], [self.shape(1), other.shape(2)]);
                end
                I = reshape(1:prod(self.siz), [self.siz(1), prod(self.siz(2:end))]);
                blktens = blktens.subsref(struct('type', {'()'}, 'subs', {{I}}));
                blktens.siz = self.siz;
                return
            end
            if ~iscell(T)
                T = {T};
            end
            blktens = other.tmprod(T, 1:length(T));
        end

        function blktens = tmprod(self, U, mode)
            % Inspired by tmprod in Tensorlab.other.subsref(S)
            %
            % This code still suffers from some issues!

            % kron matrices with I to account for block shape
            for i=1:min(length(U), 2)
                U{i} = kron(U{i}, speye(self.shape(i)));
            end

            % Numerical data
            if isnumeric(self.data)
                blktens = tmprod(self.astensor, U, mode);
                size_tens = size(blktens);
                siz = [size_tens(1:2) ./ self.shape, size_tens(3:end)];
                blktens = self.cl(blktens, siz, self.shape);
                return
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
                I = reshape(permute(I, perm), size_tens(perm(1)), []);
                S = S(I);
            end

            % Cycle through the products
            for i = 1:n
                if ~isscalar(U{i})
                    size_tens(1) = size(U{i}, 1);
                end
                S = U{i} * S;
                I = reshape(1:prod(size_tens), size_tens);
                I = permute(I, [2:N, 1]);
                if i < n
                    S = S(reshape(I, size_tens(2), []));
                    size_tens = size_tens([2:N 1]);
                end
            end

            % Inverse permute
            iperm(perm([n:N 1:n-1])) = 1:N;
            I = reshape(1:prod(size_tens), size_tens);
            I = reshape(permute(I, iperm), size_tens(iperm(1)), []);
            S = S(I);
            size_tens = size_tens(iperm);
            siz = [size_tens(1:2) ./ self.shape, size_tens(3:end)];
            blktens = self.cl(S, siz, self.shape);
        end

        function blktens = trace(self)
            % Return the trace of each of the blocks
            % TODO!
            S = struct('type', {'()', '.'}, 'subs', {{':'}, 'data'});
            blktens = self.subsref(S);
            data = blktens(1:self.shape(1):end, 1);
            for i = 2:self.shape(2)
                data = data + blktens(i:self.shape(1):end, i);
            end
            blktens = self.cl(reshape(data, self.siz), self.siz, [1, 1]);

            % S = struct('type', {'()', '.'}, 'subs', {{1}, 'data'});
            % blktens = Coefficients.repmat(trace(self.subsref(S)), self.siz);
            % for i = 1:prod(self.siz)
            %     S = struct('type', {'()', '.'}, 'subs', {{i}, 'data'});
            %     blktens = blktens.subsasgn(S(1), trace(self.subsref(S)));
            % end
        end

        function blktens = spblkdiag(self)
            % Return the entries of self on a sparse block diagonal matrix
            pr = prod(self.siz);
            i = bsxfun(@plus, repmat(1:self.shape(1), pr * self.shape(2), 1), kron(0:self.shape(2):pr*self.shape(1)-1, ones(1, self.shape(2)))')';
            j = repmat(1:pr*self.shape(2) , self.shape(1), 1);
            S = struct('type', {'()', '.'}, 'subs', {{':'}, 'data'});
            data = self.subsref(S)';
            blktens = sparse(j, i, data(:));
            % TODO: remove transposition!
            % TODO: YALMIP requires third argument to be a vector
        end

        function blktens = vertcat(varargin)
            % Find out which input is of blktens class
            % This implementation is faulty for higher dimensions!
            for i = 1:length(varargin)
                if isa(varargin{i}, mfilename)
                    t = varargin{i};
                    break
                end
            end
            S = struct('type', {'()', '.'}, 'subs', {{(1:prod(t.siz))}, 'data'});
            data = [];
            len = 0;
            for i = 1:length(varargin)
                if isa(varargin{i}, mfilename)
                    % varargin{i} = varargin{i}';
                    data = [data; varargin{i}.subsref(S)];
                    len = len + varargin{i}.shape(1);
                else
                    data = [data; repmat(varargin{i}, 1, prod(t.siz))];
                    len = len + size(varargin{i}, 1);
                end
            end
            blktens = t.cl(data, [1, prod(t.siz)], [len, t.shape(2)]);
            I = reshape(1:prod(t.siz), t.siz);
            blktens = blktens.subsref(struct('type', {'()'}, 'subs', {{I}}));
            %blktens.subsref(struct('type', {'()', '.'}, 'subs', {{':'}, 'data'}))
            %blktens.siz = t.siz

            % blktens = t.cl(zeros(size(data)), t.siz, [len, t.shape(2)]);
            % S = struct('type', {'()'}, 'subs', {{':'}});
            % blktens = blktens.subsasgn(S, data);
            % data
        end

        function blktens = horzcat(varargin)
            % Find out which input is of blktens class
            %
            % TODO: fix vertcat and horzcat!
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
                    data = [data varargin{i}.subsref(S)];
                    len = len + varargin{i}.shape(2);
                else
                    data = [data repmat(varargin{i}, prod(t.siz), 1)];
                    len = len + size(varargin{i}, 2);
                end
            end
            blktens = t.cl(data, [prod(t.siz), 1], [t.shape(1), len]);
            I = reshape(1:prod(t.siz), t.siz);
            blktens = blktens.subsref(struct('type', {'()'}, 'subs', {{I}}));
            %blktens.subsref(struct('type', {'()', '.'}, 'subs', {{':'}, 'data'}))
            %blktens.siz = t.siz
        end

        function blktens = ctranspose(self)
            % Transpose each of the blocks separately
            % blktens = self.cl(zeros(size(self.data)), self.siz, [self.shape(2), self.shape(1)]);
            blktens = self;
            blktens.shape = [self.shape(2), self.shape(1)];
            for i = 1:prod(self.siz)
                S = struct('type', {'()', '.'}, 'subs', {{i}, 'data'});
                blktens = blktens.subsasgn(S(1), self.subsref(S)');
            end
        end

        function n = end(self)
            n = prod(self.siz)
        end

        function varargout = subsref(self, S)
            % There are some serious issues here!!!!!
            switch S(1).type
                case '()'
                    if length(S(1).subs) == 1  % Linear indexing
                        if strcmp(S(1).subs{1}, ':')
                            S(1).subs{1} = (1:prod(self.siz))';
                        end
                        ssubs = size(S(1).subs{1});
                        S(1).subs{1} = reshape(S(1).subs{1}, ssubs(1), prod(ssubs(2:end)));
                        [i, j] = ind2sub([self.siz(1), prod(self.siz(2:end))], S(1).subs{1});
                        si = size(i); sj = size(j);
                        I = repmat((1:self.shape(1))', [si(1), si(2) * self.shape(2)]);
                        J = repmat(1:self.shape(2), [sj(1) * self.shape(1), sj(2)]);
                        i = kron(i, ones(self.shape));
                        j = kron(j, ones(self.shape));
                        idx = sub2ind(size(self.data), (i - 1) * self.shape(1) + I, (j - 1) * self.shape(2) + J);
                        % siz = size(S(1).subs{1});
                        y = self.cl(self.data(idx), ssubs, self.shape);
                    else % Convert indices to linear indices?
                        if length(S(1).subs) == length(self.size)
                            size_tens = self.totalsize;
                            I = reshape(1:prod(self.size), self.size);
                            I = I(S(1).subs{:});
                            I = reshape(I, size(S(1).subs{1}, 1), []);

                            %idx = sub2ind(self.siz, S(1).subs{:})
                            siz = cellfun(@length, S(1).subs);
                            S(1).subs = {I};
                            [varargout{1:nargout}] = subsref(self, S);

                            % Correct size of output
                            out = varargout{1};
                            varargout{1} = self.cl(out.data, siz, out.shape);
                            % varargout{1}.siz = siz;
                        end
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

        function c = value(self)
            c = self.cl(value(self.data), self.size, self.shape);
        end
    end

    methods (Static)
        function c = repmat(data, dims)
            % Repeats data along dims to
            shape = size(data);
            data = repmat(data, [dims(1) prod(dims(2:end))]);
            cl = str2func(mfilename);
            c = cl(data, dims, shape);
        end
    end
end
