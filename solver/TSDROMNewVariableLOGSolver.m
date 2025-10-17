classdef TSDROMNewVariableLOGSolver < TSDROMNewVariableSolverExtrapolation
    % Solver for solving a mixed-integer linear programming problem in order to generate new variables in the dual formulation of the 
    % two-stage distributionally robust optimization problem with marginal constraints. The MILP problem is formulated via the LOG
    % model.
   
    methods(Access = protected)
        function initializeModel(obj)
            % Generate the Gurobi MILP model stored in obj.Runtime.Model
            
            obj.initializeKantorovichPotentials();

            % The decision variables in the LOG formulation:
            % 1. the decision variable for the dualized second-stage problem (length = obj.problem.st2d_deci_num)
            % 2. the variables for interpolating between knots (length = size(obj.primal_solution.Kpot_duals{marg_id}, 1) + 2 for each
            %    marg_id)
            % 3. the binary variables for indicating the sub-intervals (length = ceil(log2(size( ...
            %    obj.primal_solution.Kpot_duals{marg_id}, 1) + 1)) for each marg_id)

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
                bin_num = ceil(log2(knot_num - 1));

                interp_indices{marg_id} = var_counter + (1:knot_num)';
                var_counter = var_counter + knot_num;
                bin_indices{marg_id} = var_counter + (1:bin_num)';
                var_counter = var_counter + bin_num;
                aux_lb{1, marg_id} = zeros(knot_num, 1);
                aux_lb{2, marg_id} = -inf(bin_num, 1);
                aux_ub{1, marg_id} = inf(knot_num, 1);
                aux_ub{2, marg_id} = inf(bin_num, 1);
                aux_vtype{1, marg_id} = repmat('C', knot_num, 1);
                aux_vtype{2, marg_id} = repmat('B', bin_num, 1);
            end

            model = struct;
            model.modelsense = 'max';
            model.objcon = 0;
            model.obj = zeros(var_counter, 1);
            model.obj(st2d_indices) = obj.problem.st2d_obj_act * obj.primal_solution.st1_deci + obj.problem.st2d_obj_itc;

            for marg_id = 1:obj.problem.unc_num
                model.obj(interp_indices{marg_id}) = obj.Runtime.Kantorovich_potentials{marg_id}(:, 2);
            end


            model.lb = [st2d_lb; vertcat(aux_lb{:})];
            model.ub = [st2d_ub; vertcat(aux_ub{:})];
            model.vtype = [st2d_vtype; vertcat(aux_vtype{:})];
            model.var_indices = struct;
            model.var_indices.st2d = st2d_indices;
            model.var_indices.interp = interp_indices;
            model.var_indices.bin = bin_indices;

            % The constraints in the LOG formulation:
            % 1. equality constraints requiring the sum of the interpolating variables to be equal to 1
            % 2. inequality constraints linking the binary variables and the interpolating variables
            % 3. equality constraints requiring the weighted sum of the interpolating variables to be equal to the projected dualized
            %    second-stage variable
            % 4. equality constraints characterizing the feasible set of the dualized second-stage problem

            A_Kpot_cell = cell(obj.problem.unc_num, 1);
            rhs_Kpot_cell = cell(obj.problem.unc_num, 1);
            sense_Kpot_cell = cell(obj.problem.unc_num, 1);
            constr_sum_indices = zeros(obj.problem.unc_num, 1);
            constr_counter = 0;

            for marg_id = 1:obj.problem.unc_num
                knot_num = size(obj.Runtime.Kantorovich_potentials{marg_id}, 1);
                bisect_cell = obj.compute1DBisection(knot_num - 1);
                bin_num = size(bisect_cell, 1);

                A_eq_normalize = sparse( ...
                    ones(knot_num, 1), ...
                    interp_indices{marg_id}, ...
                    1, ...
                    1, var_counter);
                rhs_eq_normalize = 1;
                sense_eq_normalize = '=';
                constr_counter = constr_counter + 1;

                A_ineq_link_cell = cell(bin_num, 1);
                rhs_ineq_link_cell = cell(bin_num, 1);
                sense_ineq_link_cell = cell(bin_num, 1);

                for bit_id = 1:bin_num
                    bit0_num = length(bisect_cell{bit_id, 1});
                    bit1_num = length(bisect_cell{bit_id, 2});

                    A_ineq_link0 = sparse( ...
                        ones(bit0_num + 1, 1), ...
                        [interp_indices{marg_id}(bisect_cell{bit_id, 1}); bin_indices{marg_id}(bit_id)], ...
                        [ones(bit0_num, 1); -1], ...
                        1, var_counter);
                    rhs_ineq_link0 = 0;
                    sense_ineq_link0 = '<';
                    constr_counter = constr_counter + 1;

                    A_ineq_link1 = sparse( ...
                        ones(bit1_num + 1, 1), ...
                        [interp_indices{marg_id}(bisect_cell{bit_id, 2}); bin_indices{marg_id}(bit_id)], ...
                        ones(bit1_num + 1, 1), ...
                        1, var_counter);
                    rhs_ineq_link1 = 1;
                    sense_ineq_link1 = '<';
                    constr_counter = constr_counter + 1;

                    A_ineq_link_cell{bit_id} = [A_ineq_link0; A_ineq_link1];
                    rhs_ineq_link_cell{bit_id} = [rhs_ineq_link0; rhs_ineq_link1];
                    sense_ineq_link_cell{bit_id} = [sense_ineq_link0; sense_ineq_link1];
                end

                A_eq_sum = sparse( ...
                    ones(obj.problem.st2d_deci_num + knot_num, 1), ...
                    [st2d_indices; interp_indices{marg_id}], ...
                    [obj.problem.st2d_obj_unc(:, marg_id); -obj.Runtime.Kantorovich_potentials{marg_id}(:, 1)], ...
                    1, var_counter);
                rhs_eq_sum = 0;
                sense_eq_sum = '=';
                constr_sum_indices(marg_id) = constr_counter + 1;
                constr_counter = constr_counter + 1;

                A_Kpot_cell{marg_id} = [A_eq_normalize; vertcat(A_ineq_link_cell{:}); A_eq_sum];
                rhs_Kpot_cell{marg_id} = [rhs_eq_normalize; vertcat(rhs_ineq_link_cell{:}); rhs_eq_sum];
                sense_Kpot_cell{marg_id} = [sense_eq_normalize; vertcat(sense_ineq_link_cell{:}); sense_eq_sum];
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

        function bisect_cell = compute1DBisection(obj, intv_num)
            % Compute a bisection of knots in a one-dimensional continuous piece-wise affine (CPWA) function used for formulating it
            % into a mixed-integer programming problem
            % Input:
            %   intv_num: number of intervals in the CPWA function
            % Output:
            %   bisect_cell: cell array containing the bisection, the number of rows is equal to ceil(log2(intv_num)) and the number of
            %   columns is 2

            bit_len = ceil(log2(intv_num));

            if ~isfield(obj.Storage, 'ReflexiveBinarySequences') || isempty(obj.Storage.ReflexiveBinarySequences)
                obj.Storage.ReflexiveBinarySequences = cell(20, 1);
            end

            if bit_len == 0
                bisect_cell = cell(0, 2);
                return;
            end

            if bit_len > 20
                error('the number of intervals is too large');
            end

            if isempty(obj.Storage.ReflexiveBinarySequences{bit_len})
                seq_mat = zeros(2^bit_len, bit_len);
                
                % fill the first two rows
                seq_mat(1:2, end) = [0; 1];

                for bit_id = 2:bit_len
                    % flip the previous matrix upside down and append to the end while adding a column of 1's to the left
                    prev_row_num = 2 ^ (bit_id - 1);
                    seq_mat(prev_row_num + (1:prev_row_num), ...
                        end - bit_id + 1:end) ...
                        = [ones(prev_row_num, 1), flipud( ...
                        seq_mat(1:prev_row_num, end - bit_id + 2:end))];
                end

                obj.Storage.ReflexiveBinarySequences{bit_len} = seq_mat;
            end

            seq_mat = obj.Storage.ReflexiveBinarySequences{bit_len}(1:intv_num, :);

            bisect_cell = cell(bit_len, 2);

            for bit_id = 1:bit_len
                rflx_col = seq_mat(:, bit_id);
                bisect_cell{bit_id, 1} = find([rflx_col(1) == 0; rflx_col(1:end - 1) == 0 & rflx_col(2:end) == 0; rflx_col(end) == 0]);
                bisect_cell{bit_id, 2} = find([rflx_col(1) == 1; rflx_col(1:end - 1) == 1 & rflx_col(2:end) == 1; rflx_col(end) == 1]);
            end
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

                if interp(1) > 1 - obj.Options.extrapolation_boundary_tolerance
                    status(marg_id, 1) = true;
                end

                if interp(end) > 1 - obj.Options.extrapolation_boundary_tolerance
                    status(marg_id, 2) = true;
                end
            end
        end
    end
end
