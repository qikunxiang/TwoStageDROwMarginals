classdef TSDROMUncertain < handle
    % Class representing the uncertain parameters in the two-stage distributionally robust optimization (DRO) problem with marginal 
    % constraints
    
    properties(SetAccess = protected, GetAccess = public)
        % Cell array containing the marginal constraints
        marginals;

        % Number of uncertain parameters
        num;

        % Matrix with two columns containing the lower and upper bounds of the supports of the marginals
        marginal_bounds = [];

        % Problem instance
        problem = [];

    end
    
    methods(Access = public)
        function obj = TSDROMUncertain(marginals)
            % Constructor
            obj.marginals = marginals;
            obj.num = length(marginals);

            obj.initializeMarginalBounds();
        end

        function setProblem(obj, problem)
            % Set the problem instance
            % Input:
            %   problem: object to TSDROMInstance

            obj.problem = problem;
        end

        function val = checkDualFeasibility(obj, points)
            % Check whether the given finite subcollection yields a feasible restricted dual problem
            % Input:
            %   points: matrix containing a finite subcollection of candidate points as rows

            current_dual_point_num = size(points, 1);
            deci_num = current_dual_point_num + obj.problem.st1_constr_num_in + obj.problem.st1_constr_num_eq;
            dual_model = struct;
            dual_model.modelsense = 'min';

            % since we are only checking for feasibility, the objective vector can be set to 0
            dual_model.objcon = 0;
            dual_model.obj = zeros(deci_num, 1);

            dual_model.lb = [zeros(current_dual_point_num, 1); -inf(obj.problem.st1_constr_num_in + obj.problem.st1_constr_num_eq, 1)];
            dual_model.ub = [inf(current_dual_point_num, 1); zeros(obj.problem.st1_constr_num_in, 1); ...
                inf(obj.problem.st1_constr_num_eq, 1)];

            % inequality constraint about convex combination of the dual points
            A_in = sparse([obj.problem.st2d_obj_act' * points', -obj.problem.st1_constr_mat_in', -obj.problem.st1_constr_mat_eq']);
            rhs_in = -obj.problem.st1_obj;

            % second equality constraint requiring all weights to sum to 1
            A_eq = sparse(ones(current_dual_point_num, 1), (1:current_dual_point_num)', 1, 1, deci_num);
            rhs_eq = 1;

            dual_model.A = [A_in; A_eq];
            dual_model.rhs = [rhs_in; rhs_eq];
            dual_model.sense = [repmat('>', obj.problem.st1_deci_num, 1); '='];

            params = struct('OutputFlag', 0);

            result = gurobi(dual_model, params);

            val = strcmp(result.status, 'OPTIMAL');
        end

        function dual_point = computeOneDualPoint(obj)
            % Compute a single dual point via solving an LP problem. It is guaranteed that a single dual point in the restricted dual
            % problem is sufficient to make it feasible.
            % Output:
            %   dual_point: vector representing the computed dual point

            % the decision variables in this LP problem include:
            % 1. the candidate dual point
            % 2. the dual variable associated with the first-stage inequality constraints
            % 3. the dual variable associated with the first-stage equality constraints
            deci_num = obj.problem.st2d_deci_num + obj.problem.st1_constr_num_in + obj.problem.st1_constr_num_eq;

            model = struct;
            model.modelsense = 'min';

            % Since we only want to find a feasible point of the LP problem, the objective is set to 0
            model.objcon = 0;
            model.obj = zeros(deci_num, 1);

            model.lb = [obj.problem.st2d_constr_lb; -inf(obj.problem.st1_constr_num_in, 1); ...
                -inf(obj.problem.st1_constr_num_eq, 1)];
            model.ub = [obj.problem.st2d_constr_ub; zeros(obj.problem.st1_constr_num_in, 1); ...
                inf(obj.problem.st1_constr_num_eq, 1)];

            % the inequality constraint that relates the candidate dual point, the dual variables associated with the first-stage
            % constraints, and the first-stage objective vector
            A_in1 = sparse([obj.problem.st2d_obj_act', -obj.problem.st1_constr_mat_in', -obj.problem.st1_constr_mat_eq']);
            rhs_in1 = -obj.problem.st1_obj;
            sense_in1 = repmat('>', obj.problem.st1_deci_num, 1);

            % the inequality constraint that comes from the dualized second-stage problem
            A_in2 = [sparse(obj.problem.st2d_constr_mat_in), ...
                sparse(obj.problem.st2d_constr_num_in, obj.problem.st1_constr_num_in + obj.problem.st1_constr_num_eq)];
            rhs_in2 = obj.problem.st2d_constr_rhs_in;
            sense_in2 = repmat('<', obj.problem.st2d_constr_num_in, 1);

            model.A = [A_in1; A_in2];
            model.rhs = [rhs_in1; rhs_in2];
            model.sense = [sense_in1; sense_in2];

            params = struct('OutputFlag', 0);

            result = gurobi(model, params);

            if ~strcmp(result.status, 'OPTIMAL')
                error('the original two-stage DRO problem is unbounded from below');
            end

            dual_point = result.x(1:obj.problem.st2d_deci_num);
        end
    end

    methods(Access = protected)
        function initializeMarginalBounds(obj)
            % initialize the lower and upper bounds of the supports of the marginals

            obj.marginal_bounds = zeros(obj.num, 2);

            for marg_id = 1:obj.num
                obj.marginal_bounds(marg_id, 1) = obj.marginals{marg_id}.Supp.LowerBound;
                obj.marginal_bounds(marg_id, 2) = obj.marginals{marg_id}.Supp.UpperBound;
            end
        end
    end
end

