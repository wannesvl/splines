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
