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
            %
            % Only supports 2D splines at the moment
            %
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

        function s = f(self, x)
            s = (self.basis{2}.f(x(:, 2)) * (self.basis{1}.f(x(:, 1)) * self.coeffs)')';
            % for i = 1:self.dims()
            %     s = s;
            % end
        end

        function d = dims(self)
            d = length(self.basis)
        end

        function s = plus(self, other)
            if isa(other, 'TensorBSpline')
                basis = cellfun(@plus, self.basis, other.basis, 'UniformOutput', false)
                Tself = cellfun(@(b1, b2) b1.transform(b2), basis, self.basis, 'UniformOutput', false)
                Tother = cellfun(@(b1, b2) b1.transform(b2), basis, other.basis, 'UniformOutput', false)
                cself = self.coeffs
                for i=1:self.dims()
                    cself = permute(Tself{i} * permute(cself, 1, i), 1, i)
                end
                cother = other.coeffs
                for i=1:self.dims()
                    cother = permute(Tself{i} * permute(cother, 1, i), 1, i)
                end
                coeffs = cself + cother
            else
                basis = self.basis
                coeffs = self.coeffs + other
            end
            s = self.cl(basis, coeffs)
        end
    end
end
