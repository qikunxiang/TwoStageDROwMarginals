classdef (Abstract) TSDROMNewVariableSolverExtrapolation < TSDROMNewVariableSolver
    % Solver for solving a mixed-integer linear programming problem in order to generate new variables in the dual formulation of the 
    % two-stage distributionally robust optimization problem with marginal constraints. The solver will iteratively update the
    % extrapolation distances beyond the first and last knots in the mixed-integer linear programming formulation.
    
    methods(Access = public)

        function setProblem(obj, problem)
            % Set the problem that the solver is solving
            % Input: 
            %   problem: instance of TSDROMInstance

            setProblem@TSDROMNewVariableSolver(obj, problem);

            % initial extrapolation distance beyond the first and last knots
            if ~isfield(obj.Options, 'init_extrapolation') || isempty(obj.Options.init_extrapolation)
                obj.Options.init_extrapolation = 10;
            end

            % lower bound for the leftmost knot; once the current extrapolation reaches this lower bound, no further extrapolation is
            % performed
            if ~isfield(obj.Options, 'proj_lb') || isempty(obj.Options.proj_lb)
                obj.Options.proj_lb = -5e3 * ones(obj.problem.unc_num, 1);
            end

            % upper bound for the rightmost knot; once the current extrapolation reaches this upper bound, no further extrapolation is
            % performed
            if ~isfield(obj.Options, 'proj_ub') || isempty(obj.Options.proj_ub)
                obj.Options.proj_ub = 5e3 * ones(obj.problem.unc_num, 1);
            end

            % tolerance value as ratio for determining whether a solution is on the boundary of extrapolation
            if ~isfield(obj.Options, 'extrapolation_boundary_tolerance') ...
                    || isempty(obj.extrapolation_boundary_tolerance)
                obj.Options.extrapolation_boundary_tolerance = 1e-3;
            end

            % tolerance value for determining whether a solution is on the boundary of extrapolation
            if ~isfield(obj.Options, 'boundary_tolerance') || isempty(obj.boundary_tolerance)
                obj.Options.boundary_tolerance = 1e-5;
            end

            obj.initializeExtrapolation();
        end
        
        function [pool_st2d, gap, output] = run(obj)
            % Run the solver
            % Outputs:
            %   pool_st2d: matrix containing the computed approximate minimizers in the feasible set of the dualized second-stage
            %   problem
            %   gap: computed sub-optimality gap
            %   output: struct containing additional output information, such as trial_num, total_time

            % start the timer
            total_timer = tic;

            if isempty(obj.problem)
                error('problem must be set first');
            end

            if isempty(obj.primal_solution)
                error('the primal solution must be set first');
            end

            % open the log file
            if ~isempty(obj.Options.log_file)
                log_file = fopen(obj.Options.log_file, 'a');

                if log_file < 0
                    error('cannot open log file');
                end

                fprintf(log_file, '--- new variable solver starts ---\n');
            end

            % options of the mixed-integer solver
            gl_options = struct;
            gl_options.OutputFlag = 0;
            gl_options.IntFeasTol = 1e-6;
            gl_options.FeasibilityTol = 1e-8;
            gl_options.OptimalityTol = 1e-8;
            gl_options.PoolSolutions = 100;
            gl_options.PoolGap = 0.8;
            gl_options.NodefileStart = 2;
            gl_options.BestBdStop = -1e-6;
            gl_options.BestObjStop = inf;
            gl_options.MIPGap = 1e-4;
            gl_options.MIPGapAbs = 1e-10;

            % set the additional options for the mixed-integer solver
            gl_options_fields = fieldnames(obj.GlobalOptions);
            gl_options_values = struct2cell(obj.GlobalOptions);

            for fid = 1:length(gl_options_fields)
                gl_options.(gl_options_fields{fid}) = gl_options_values{fid};
            end

            trial_num = 0;

            while true
                trial_num = trial_num + 1;
                model = obj.Runtime.Model;
    
                result = gurobi(model, gl_options);
    
                if ~strcmp(result.status, 'OPTIMAL') && ~strcmp(result.status, 'USER_OBJ_LIMIT') ...
                        && ~strcmp(result.status, 'TIME_LIMIT')
                    error('error in the mixed-integer solver');
                end

                boundary_status = obj.checkBoundaryStatus(result.x);

                if all(all(~boundary_status))
                    gap = result.objbound;
        
                    % get a set of approximate optimizers
                    if isfield(result, 'pool')
                        pool_sols = horzcat(result.pool.xn)';
                    else
                        pool_sols = result.x';
                    end
        
                    pool_st2d = pool_sols(:, model.var_indices.st2d);
                    pool_st2d = obj.processNewVariables(pool_st2d);

                    % display
                    if obj.Options.display
                        fprintf('successfully generated %d new variables, gap = %.6f\n', size(pool_st2d, 1), gap);
                    end

                    % logging
                    if ~isempty(obj.Options.log_file)
                        fprintf(log_file, 'successfully generated %d new variables, gap = %.6f\n', size(pool_st2d, 1), gap);
                    end

                    break;
                else
                    % display
                    if obj.Options.display
                        fprintf('unsuccessful due to extrapolation range being too narrow\n');
                    end

                    % logging
                    if ~isempty(obj.Options.log_file)
                        fprintf(log_file, 'unsuccessful due to extrapolation range being too narrow\n');
                    end

                    for marg_id = 1:obj.problem.unc_num
                        if boundary_status(marg_id, 1)
                            new_extrapolation_left = obj.Runtime.Extrapolations(marg_id, 1) * 2;

                            % display
                            if obj.Options.display
                                fprintf('    %5s extrapolation distance of dimension %3d increased from %.0f to %.0f\n', ...
                                    'left', marg_id, obj.Runtime.Extrapolations(marg_id, 1), new_extrapolation_left);
                            end

                            % logging
                            if ~isempty(obj.Options.log_file)
                                fprintf(log_file, '    %5s extrapolation distance of dimension %3d increased from %.0f to %.0f\n', ...
                                    'left', marg_id, obj.Runtime.Extrapolations(marg_id, 1), new_extrapolation_left);
                            end

                            obj.Runtime.Extrapolations(marg_id, 1) = new_extrapolation_left;
                        end

                        if boundary_status(marg_id, 2)
                            new_extrapolation_right = obj.Runtime.Extrapolations(marg_id, 2) * 2;
                            
                            % display
                            if obj.Options.display
                                fprintf('    %5s extrapolation distance of dimension %3d increased from %f to %f\n', ...
                                    'right', marg_id, obj.Runtime.Extrapolations(marg_id, 2), new_extrapolation_right);
                            end

                            % logging
                            if ~isempty(obj.Options.log_file)
                                fprintf(log_file, '    %5s extrapolation distance of dimension %3d increased from %f to %f\n', ...
                                    'right', marg_id, obj.Runtime.Extrapolations(marg_id, 2), new_extrapolation_right);
                            end

                            obj.Runtime.Extrapolations(marg_id, 2) = new_extrapolation_right;
                        end
                    end

                    obj.initializeModel();
                end
            end

            % close the log file
            if ~isempty(obj.Options.log_file)
                fprintf(log_file, '--- new variable solver ends ---\n\n');
                fclose(log_file);
            end

            output = struct;
            output.trial_num = trial_num;

            % stop the timer
            total_time = toc(total_timer);
            output.total_time = total_time;
        end
    end

    methods(Access = protected)

        function initializeExtrapolation(obj)
            % Initialize the extrapolation distances beyond the first and the last knots, as well as the lower/upper bounds
            obj.Runtime.Extrapolations = obj.Options.init_extrapolation * ones(obj.problem.unc_num, 2);

            obj.Runtime.ProjLowerBounds = max(obj.problem.st2d_proj_lb, obj.Options.proj_lb);
            obj.Runtime.ProjUpperBounds = min(obj.problem.st2d_proj_ub, obj.Options.proj_ub);
        end

        function initializeKantorovichPotentials(obj)
            % Computes the knots and function values in the Kantorovich potential functions

            obj.Runtime.Kantorovich_potentials = cell(obj.problem.unc_num, 1);

            for marg_id = 1:obj.problem.unc_num
                Kpot_mat = obj.primal_solution.Kpot_duals_cell{marg_id};
                knots = Kpot_mat(:, 1);
                vals = Kpot_mat(:, 2);

                % if the leftmost knot is not equal to the given lower bound
                if knots(1) > obj.Runtime.ProjLowerBounds(marg_id) + obj.Options.boundary_tolerance
                    new_knot_left = max(knots(1) - obj.Runtime.Extrapolations(marg_id, 1), obj.problem.st2d_proj_lb(marg_id));
                    new_val_left = vals(1) - (knots(1) - new_knot_left) * Kpot_mat(1, 3);

                    knots = [new_knot_left; knots]; %#ok<AGROW>
                    vals = [new_val_left; vals]; %#ok<AGROW>
                end

                % if the rightmost knot is not equal to the given upper bound
                if knots(end) < obj.Runtime.ProjUpperBounds(marg_id) - obj.Options.boundary_tolerance
                    new_knot_right = min(knots(end) + obj.Runtime.Extrapolations(marg_id, 2), obj.problem.st2d_proj_ub(marg_id));
                    new_val_right = vals(end) + (new_knot_right - knots(end)) * Kpot_mat(end, 4);

                    knots = [knots; new_knot_right]; %#ok<AGROW>
                    vals = [vals; new_val_right]; %#ok<AGROW>
                end

                obj.Runtime.Kantorovich_potentials{marg_id} = [knots, vals];
            end
        end

        function status = checkBoundaryStatus(obj, sol)
            % Check whether the extrapolation distance at each boundary needs to be doubled
            % Input:
            %   sol: column vector containing the solution to be checked
            % Output:
            %   status: two-column logical matrix indicating whether the left or right boundary of extrapolation for each marginal
            %   needs to be doubled

            status = obj.checkExtrapolationBounds(sol);

            for marg_id = 1:obj.problem.unc_num
                % if the current leftmost/rightmost knot is already at the lower/upper bound, do not increase extrapolation
                if obj.Runtime.Kantorovich_potentials{marg_id}(1, 1) <= obj.Runtime.ProjLowerBounds(marg_id) ...
                        + obj.Options.boundary_tolerance
                    status(marg_id, 1) = false;
                end

                if obj.Runtime.Kantorovich_potentials{marg_id}(end, 1) >= obj.Runtime.ProjUpperBounds(marg_id) ...
                        - obj.Options.boundary_tolerance
                    status(marg_id, 2) = false;
                end
            end
        end
    end

    methods(Abstract, Access = protected)

        % Check from a solution of the MILP problem whether it is on the boundary of extrapolation
        status = checkExtrapolationBounds(obj, sol);
    end
end

