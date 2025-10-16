classdef TSDROMNewVariableINCSolver < TSDROMNewVariableSolverExtrapolation
    % Solver for solving a mixed-integer linear programming problem in order to generate new variables in the dual formulation of the 
    % two-stage distributionally robust optimization problem with marginal constraints. The MILP problem is formulated via the INC
    % model.

    methods(Access = protected)
        function initializeModel(obj)
            % Generate the Gurobi MILP model stored in obj.Runtime.Model

            obj.initializeKantorovichPotentials();

            % The decision variables in the INC formulation:
            % 1. the decision variable for the dualized second-stage problem (length = obj.problem.st2d_deci_num)
            % 2. the variables for interpolating between knots (length = size(obj.primal_solution.Kpot_duals_cell{marg_id}, 1) + 1 for 
            %    each marg_id)
            % 3. the binary variables for indicating the sub-intervals (length = size(obj.primal_solution.Kpot_duals_cell{marg_id}, 1) 
            %    for each marg_id)

            var_counter = 0;
            st2d_indices = var_counter + (1:obj.problem.st2d_deci_num)';
            var_counter = var_counter + obj.problem.st2d_deci_num;
            st2d_lb = obj.problem.st2d_constr_lb;
            st2d_ub = obj.problem.st2d_constr_ub;
            st2d_vtype = repmat('C', length(st2d_indices), 1);

            interp_indices = cell(obj.problem.unc_num, 1);
            bin_indices = cell(obj.problem.unc_num, 1);
            aux_lb = cell(2, obj.problem.unc_num);
            aux_ub = cell(2, obj.problem.unc_num);
            aux_vtype = cell(2, obj.problem.unc_num);

            for marg_id = 1:obj.problem.unc_num
                knot_num = size(obj.Runtime.Kantorovich_potentials{marg_id}, 1);

                interp_indices{marg_id} = var_counter + (1:knot_num - 1)';
                var_counter = var_counter + knot_num - 1;
                bin_indices{marg_id} = var_counter + (1:knot_num - 2)';
                var_counter = var_counter + knot_num - 2;
                aux_lb{1, marg_id} = [-inf(knot_num - 2, 1); 0];
                aux_lb{2, marg_id} = -inf(knot_num - 2, 1);
                aux_ub{1, marg_id} = [1; inf(knot_num - 2, 1)];
                aux_ub{2, marg_id} = inf(knot_num - 2, 1);
                aux_vtype{1, marg_id} = repmat('C', knot_num - 1, 1);
                aux_vtype{2, marg_id} = repmat('B', knot_num - 2, 1);
            end

            model = struct;
            model.modelsense = 'max';
            model.objcon = 0;
            model.obj = zeros(var_counter, 1);
            model.obj(st2d_indices) = obj.problem.st2d_obj_act * obj.primal_solution.st1_deci + obj.problem.st2d_obj_itc;

            for marg_id = 1:obj.problem.unc_num
                model.objcon = model.objcon + obj.Runtime.Kantorovich_potentials{marg_id}(1, 2);
                model.obj(interp_indices{marg_id}) = diff(obj.Runtime.Kantorovich_potentials{marg_id}(:, 2));
            end
            
            model.lb = [st2d_lb; vertcat(aux_lb{:})];
            model.ub = [st2d_ub; vertcat(aux_ub{:})];
            model.vtype = [st2d_vtype; vertcat(aux_vtype{:})];
            model.var_indices = struct;
            model.var_indices.st2d = st2d_indices;
            model.var_indices.interp = interp_indices;
            model.var_indices.bin = bin_indices;

            % The constraints in the INC formulation:
            % 1. inequality constraints linking the binary variables and the interpolating variables
            % 2. equality constraints requiring the weighted sum of the interpolating variables to be equal to the projected dualized
            %    second-stage variable
            % 3. inequality constraints characterizing the feasible set of the dualized second-stage problem

            A_Kpot_cell = cell(obj.problem.unc_num, 1);
            rhs_Kpot_cell = cell(obj.problem.unc_num, 1);
            sense_Kpot_cell = cell(obj.problem.unc_num, 1);
            constr_sum_indices = zeros(obj.problem.unc_num, 1);
            constr_counter = 0;

            for marg_id = 1:obj.problem.unc_num
                knot_num = size(obj.Runtime.Kantorovich_potentials{marg_id}, 1);

                A_ineq_link1 = sparse( ...
                    repmat((1:knot_num - 2)', 2, 1), ...
                    [interp_indices{marg_id}(1:end - 1); bin_indices{marg_id}], ...
                    [ones(knot_num - 2, 1); -ones(knot_num - 2, 1)], ...
                    knot_num - 2, var_counter);
                rhs_ineq_link1 = zeros(knot_num - 2, 1);
                sense_ineq_link1 = repmat('>', knot_num - 2, 1);
                constr_counter = constr_counter + knot_num - 2;

                A_ineq_link2 = sparse( ...
                    repmat((1:knot_num - 2)', 2, 1), ...
                    [interp_indices{marg_id}(2:end); bin_indices{marg_id}], ...
                    [ones(knot_num - 2, 1); -ones(knot_num - 2, 1)], ...
                    knot_num - 2, var_counter);
                rhs_ineq_link2 = zeros(knot_num - 2, 1);
                sense_ineq_link2 = repmat('<', knot_num - 2, 1);
                constr_counter = constr_counter + knot_num - 2;

                A_eq_sum = sparse( ...
                    ones(obj.problem.st2d_deci_num + knot_num - 1, 1), ...
                    [st2d_indices; interp_indices{marg_id}], ...
                    [obj.problem.st2d_obj_unc(:, marg_id); -diff(obj.Runtime.Kantorovich_potentials{marg_id}(:, 1))], ...
                    1, var_counter);
                rhs_eq_sum = obj.Runtime.Kantorovich_potentials{marg_id}(1, 1);
                sense_eq_sum = '=';
                constr_sum_indices(marg_id) = constr_counter + 1;
                constr_counter = constr_counter + 1;

                A_Kpot_cell{marg_id} = [A_ineq_link1; A_ineq_link2; A_eq_sum];
                rhs_Kpot_cell{marg_id} = [rhs_ineq_link1; rhs_ineq_link2; rhs_eq_sum];
                sense_Kpot_cell{marg_id} = [sense_ineq_link1; sense_ineq_link2; sense_eq_sum];
            end

            A_st2d = sparse( ...
                repmat((1:obj.problem.st2d_constr_num_in)', obj.problem.st2d_deci_num, 1), ...
                repelem(st2d_indices, obj.problem.st2d_constr_num_in, 1), ...
                obj.problem.st2d_constr_mat_in(:), ...
                obj.problem.st2d_constr_num_in, var_counter);
            rhs_st2d = obj.problem.st2d_constr_rhs_in;
            sense_st2d = repmat('<', obj.problem.st2d_constr_num_in, 1);

            model.A = [vertcat(A_Kpot_cell{:}); A_st2d];
            model.rhs = [vertcat(rhs_Kpot_cell{:}); rhs_st2d];
            model.sense = [vertcat(sense_Kpot_cell{:}); sense_st2d];
            model.constr_sum_indices = constr_sum_indices;

            obj.Runtime.Model = model;
        end

        function status = checkExtrapolationBounds(obj, sol)
            % Check from a solution of the MILP problem whether it is on the boundary of extrapolation
            % Input:
            %   sol: column vector containing the solution to be checked
            % Output:
            %   status: two-column logical matrix indicating whether the solution is on the left or right boundary of extrapolation for
            %   each marginal

            status = false(obj.problem.unc_num, 2);
            var_indices = obj.Runtime.Model.var_indices;

            for marg_id = 1:obj.problem.unc_num
                interp = sol(var_indices.interp{marg_id});

                if interp(1) < obj.Options.extrapolation_boundary_tolerance
                    status(marg_id, 1) = true;
                end

                if interp(end) > 1 - obj.Options.extrapolation_boundary_tolerance
                    status(marg_id, 2) = true;
                end
            end
        end
    end
end

