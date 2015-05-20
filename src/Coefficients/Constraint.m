
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

classdef Constraint
    properties
        expr
        lb
        ub
    end
    properties (Access=protected)
        cl
    end
    methods
        function c = Constraint(expr, lb, ub)
            c.expr = expr;
            if isscalar(lb)
                lb = lb * ones(size(expr));
            end
            if isscalar(ub)
                ub = ub * ones(size(expr));
            end
            c.lb = lb;
            c.ub = ub;
            c.cl = str2func(mfilename);
        end

        function c = vertcat(varargin)
            expr = varargin{1}.expr;
            lb = varargin{1}.lb;
            ub = varargin{1}.ub;
            for i = 2:length(varargin)
                expr = vertcat(expr, varargin{i}.expr);
                lb = [lb; varargin{i}.lb];
                ub = [ub; varargin{i}.ub];
            end
            c = varargin{1}.cl(expr, lb, ub);
        end

        function c = horzcat(varargin)
            c = vertcat(varargin{:});
        end
    end
end
