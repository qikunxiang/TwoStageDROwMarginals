classdef TSDROMNewVariablePWLSolver < TSDROMNewVariableSolver
    % Solver for solving a mixed-integer linear programming problem in order to generate new variables in the dual formulation of the 
    % two-stage distributionally robust optimization problem with marginal constraints. The problem is formulated directly via Gurobi's
    % piece-wise linear functionalities

    methods(Access = public)
        function obj = TSDROMNewVariablePWLSolver(varargin)
            % Constructor
            obj@TSDROMNewVariableSolver(varargin{:});

            % set the default option for the maximum discrepancy in the left and right derivatives of a Kantorovich potential function
            % for it to be considered affine at a knot
            if ~isfield(obj.Options, 'Kpot_slope_tolerance') || isempty(obj.Options.Kpot_slope_tolerance)
                obj.Options.Kpot_slope_tolerance = 1e-6;
            end

            % set the default option for the tolerance value for a knot to be regarded as being on the lower/upper boundary
            if ~isfield(obj.Options, 'boundary_tolerance') || isempty(obj.Options.boundary_tolerance)
                obj.Options.boundary_tolerance = 1e-3;
            end
        end
    end

    methods(Access = protected)
        
        function initializeKantorovichPotentials(obj)
            % Computes the knots and function values in the Kantorovich potential functions

            obj.Runtime.Kantorovich_potentials = cell(obj.problem.unc_num, 1);

            for marg_id = 1:obj.problem.unc_num
                Kpot_mat = obj.primal_solution.Kpot_duals_cell{marg_id};
                knots = Kpot_mat(:, 1);
                vals = Kpot_mat(:, 2);

                % if the leftmost knot is already equal to the lower bound, there is no need for an additional knot
                if knots(1) > obj.problem.st2d_proj_lb(marg_id) + obj.Options.boundary_tolerance
                    % if the left-derivative at the leftmost knot is not equal to the right-derivative, an additional knot needs to be
                    % added to encode the left ray; otherwise we move this knot closer to the next knot to the right
                    if Kpot_mat(1, 4) - Kpot_mat(1, 3) > obj.Options.Kpot_slope_tolerance
                        new_knot_left = knots(1) - 1;
                        new_val_left = vals(1) - Kpot_mat(1, 3);
                        knots = [new_knot_left; knots]; %#ok<AGROW>
                        vals = [new_val_left; vals]; %#ok<AGROW>
                    else
                        knots(1) = knots(2) - 1;
                        vals(1) = vals(2) - Kpot_mat(2, 3);
                    end
                end

                % if the rightmost knot is already equal to the upper bound, there is no need for an additional knot
                if knots(end) < obj.problem.st2d_proj_ub(marg_id) - obj.Options.boundary_tolerance
                    % if the right-derivative at the rightmost knot is not equal to the left-derivative, an additional knot needs to be
                    % added to encode the right ray; otherwise we move this knot closer to the next knot to the left
                    if Kpot_mat(end, 4) - Kpot_mat(end, 3) > obj.Options.Kpot_slope_tolerance
                        new_knot_right = knots(end) + 1;
                        new_val_right = vals(end) + Kpot_mat(end, 4);
                        knots = [knots; new_knot_right]; %#ok<AGROW>
                        vals = [vals; new_val_right]; %#ok<AGROW>
                    else
                        knots(end) = knots(end - 1) + 1;
                        vals(end) = vals(end - 1) + Kpot_mat(end - 1, 4);
                    end
                end

                obj.Runtime.Kantorovich_potentials{marg_id} = [knots, vals];
            end
        end

        function initializeModel(obj)
            % Generate the Gurobi model stored in obj.Runtime.Model

            obj.initializeKantorovichPotentials();

            % The decision variables in the INC formulation:
            % 1. the decision variable for the dualized second-stage problem (length = obj.problem.st2d_deci_num)
            % 2. the variables representing the projections of the decision variable for the dualized second-stage problem (length =
            % obj.problem.unc_num)
            % 3. the variables representing the values of the piece-wise linear functions (length = obj.problem.unc_num)

            var_counter = 0;
            st2d_indices = var_counter + (1:obj.problem.st2d_deci_num)';
            var_counter = var_counter + obj.problem.st2d_deci_num;
            st2d_lb = obj.problem.st2d_constr_lb;
            st2d_ub = obj.problem.st2d_constr_ub;
            st2d_vtype = repmat('C', obj.problem.st2d_deci_num, 1);

            proj_indices = var_counter + (1:obj.problem.unc_num)';
            var_counter = var_counter + obj.problem.unc_num;
            proj_lb = obj.problem.st2d_proj_lb;
            proj_ub = obj.problem.st2d_proj_ub;
            proj_vtype = repmat('C', obj.problem.unc_num, 1);

            pwl_indices = var_counter + (1:obj.problem.unc_num)';
            var_counter = var_counter + obj.problem.unc_num;
            pwl_lb = -inf(obj.problem.unc_num, 1);
            pwl_ub = inf(obj.problem.unc_num, 1);
            pwl_vtype = repmat('C', obj.problem.unc_num, 1);

            model = struct;
            model.modelsense = 'max';
            model.objcon = 0;
            model.obj = zeros(var_counter, 1);
            model.obj(st2d_indices) = obj.problem.st2d_obj_act * obj.primal_solution.st1_deci + obj.problem.st2d_obj_itc;
            model.obj(pwl_indices) = 1;
            
            model.lb = [st2d_lb; proj_lb; pwl_lb];
            model.ub = [st2d_ub; proj_ub; pwl_ub];
            model.vtype = [st2d_vtype; proj_vtype; pwl_vtype];
            model.var_indices = struct;
            model.var_indices.st2d = st2d_indices;
            model.var_indices.proj = proj_indices;
            model.var_indices.pwl = pwl_indices;

            % The constraints in the global optimization problem:
            % 1. piece-wise linear equality constraints
            % 2. equality constraints linking the dualized second-stage decision variables and the projections
            % 3. inequality constraints characterizing the feasible set of the dualized second-stage problem

            xvar_cell = cell(obj.problem.unc_num, 1);
            yvar_cell = cell(obj.problem.unc_num, 1);
            xpts_cell = cell(obj.problem.unc_num, 1);
            ypts_cell = cell(obj.problem.unc_num, 1);

            for marg_id = 1:obj.problem.unc_num
                xvar_cell{marg_id} = proj_indices(marg_id);
                yvar_cell{marg_id} = pwl_indices(marg_id);
                xpts_cell{marg_id} = obj.Runtime.Kantorovich_potentials{marg_id}(:, 1);
                ypts_cell{marg_id} = obj.Runtime.Kantorovich_potentials{marg_id}(:, 2);
            end

            st2d_obj_unc_t = obj.problem.st2d_obj_unc';
            A_link = sparse( ...
                repmat((1:obj.problem.unc_num)', obj.problem.st2d_deci_num + 1, 1), ...
                [repelem(st2d_indices, obj.problem.unc_num, 1); proj_indices], ...
                [st2d_obj_unc_t(:); -ones(obj.problem.unc_num, 1)], ...
                obj.problem.unc_num, var_counter);
            rhs_link = zeros(obj.problem.unc_num, 1);
            sense_link = repmat('=', obj.problem.unc_num, 1);

            A_st2d = sparse( ...
                repmat((1:obj.problem.st2d_constr_num_in)', obj.problem.st2d_deci_num, 1), ...
                repelem(st2d_indices, obj.problem.st2d_constr_num_in, 1), ...
                obj.problem.st2d_constr_mat_in(:), ...
                obj.problem.st2d_constr_num_in, var_counter);
            rhs_st2d = obj.problem.st2d_constr_rhs_in;
            sense_st2d = repmat('<', obj.problem.st2d_constr_num_in, 1);

            model.genconpwl = struct('xvar', xvar_cell, 'yvar', yvar_cell, 'xpts', xpts_cell, 'ypts', ypts_cell);
            model.A = [A_link; A_st2d];
            model.rhs = [rhs_link; rhs_st2d];
            model.sense = [sense_link; sense_st2d];

            obj.Runtime.Model = model;
        end
    end
end

