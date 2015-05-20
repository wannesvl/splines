classdef TensorBasis
    properties
        basis  % A cell array of bases
    end
    properties (Access=protected)
        cl
    end

    methods
        function basis = TensorBasis(varargin)
            % Constructor for TensorBasis object
            self.basis = varargin;
            self.cl = mfilename;
        end

        function b = f(self, varargin)
            % Evaluate basis
            b = cellfun(@(b, x) b.f(x), self.basis, varargin, 'UniformOutput', false);
        end

        function basis = plus(self, other)
            b = cellfun(@(b, x) b.f(x), self.basis, other.basis, 'UniformOutput', false);
            basis = self.cl(b);
        end

        function
    end
end
