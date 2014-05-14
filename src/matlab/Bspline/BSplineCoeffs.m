classdef BSplineCoeffs
    properties (Access=protected)
        cl
        i
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
            if strcmp(class(coeffs), 'cell')
                % Check if all elements are of equal size
                sizes = cellfun(@size, coeffs, 'UniformOutpu', false);
                if isequal(sizes{:})
                    c.coeffs = coeffs;
                else
                    error('Coefficients should all be of equal size')
                end
            elseif isvector(coeffs)  % Scalar valued
                coeffs = coeffs(:);
                c.coeffs = mat2cell(coeffs, ones(size(coeffs, 1), 1));
            else
                error('Coeffs not correctly formatted')
            end
            c.i = eye(size(c.coeffs{1}, 1));
            c.cl = str2func(class(c));
        end

        function c = mtimes(a, self)
            % Implement 'generalized' inner product
            if strcmp(class(a), 'double')  % 
                c = cell2mat(self.coeffs);
                if ~isscalar(self.i)
                    c = c';
                end
                c = kron(a, self.i) * c;
            end
        end

        function c = transform(self, T)
            % Transform coefficients with transformation matrix T
            coeffs = T * self
            coeffs = mat2cell(coeffs)
            c = self.cl(coeffs)
        end

        % function c = subsref(self, s)
        %     c = self.coeffs{s}
        % end
    end
end