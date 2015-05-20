
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
