classdef (Abstract) TSDROMNewVariableSolver < handle
    % Solver for solving a mixed-integer linear programming problem in order to generate new variables in the dual formulation of the 
    % two-stage distributionally robust optimization problem with marginal constraints
    
    properties(SetAccess = protected, GetAccess = public)
        Options = [];

        GlobalOptions = [];

        problem = [];

        primal_solution = [];

        Storage = struct;
        Runtime = struct;
    end
    
    methods(Access = public)
        function obj = TSDROMNewVariableSolver(opt, gl_opt)
            % Constructor
            % Inputs: 
            %   problem: instance of TSDROMInstance
            %   opt: struct containing the options for the new variable solver
            %   gl_opt: struct containing the options for the Gurobi MIP solver
            
            if ~exist('opt', 'var') || isempty(opt)
                opt = struct;
            end

            obj.Options = opt;

            % display
            if ~isfield(obj.Options, 'display') || isempty(obj.Options.display)
                obj.Options.display = false;
            end

            % log file
            if ~isfield(obj.Options, 'log_file') || isempty(obj.Options.log_file)
                obj.Options.log_file = '';
            end

            if ~exist('gl_opt', 'var') || isempty(gl_opt)
                gl_opt = struct;
            end

            obj.GlobalOptions = gl_opt;
        end

        function updateOptions(obj, opt, gl_opt)
            % Update the options for the solver and the global MILP solver
            % Inputs:
            %   opt: struct containing the additional options for the solver
            %   gl_opt: struct containing the additional options for the global MILP solver
            
            opt_fields = fieldnames(opt);
            opt_values = struct2cell(opt);

            for fid = 1:length(opt_fields)
                obj.Options.(opt_fields{fid}) = opt_values{fid};
            end

            gl_opt_fields = fieldnames(gl_opt);
            gl_opt_values = struct2cell(gl_opt);

            for fid = 1:length(gl_opt_fields)
                obj.GlobalOptions.(gl_opt_fields{fid}) = gl_opt_values{fid};
            end
        end

        function setProblem(obj, problem)
            % Set the problem that the solver is solving
            % Input: 
            %   problem: instance of TSDROMInstance

            if ~isempty(obj.problem)
                error('problem has already been set');
            end

            obj.problem = problem;

            obj.primal_solution = [];
            obj.Storage = struct;
            obj.Runtime = struct;
        end

        function setPrimalSolution(obj, primal_sol)
            % Set the primal solution
            % Input:
            %   primal_solution: struct representing the current primal optimizer computed by the dual LSIP solver
            
            if isempty(obj.problem)
                error('problem must be set first');
            end

            obj.primal_solution = primal_sol;

            obj.initializeModel();
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

            % set the additional options for the mixed-integer solver
            gl_options_fields = fieldnames(obj.GlobalOptions);
            gl_options_values = struct2cell(obj.GlobalOptions);

            for fid = 1:length(gl_options_fields)
                gl_options.(gl_options_fields{fid}) = gl_options_values{fid};
            end

            model = obj.Runtime.Model;

            result = gurobi(model, gl_options);

            if ~strcmp(result.status, 'OPTIMAL') && ~strcmp(result.status, 'USER_OBJ_LIMIT') ...
                        && ~strcmp(result.status, 'TIME_LIMIT')
                error('error in the mixed-integer solver');
            end


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

            % close the log file
            if ~isempty(obj.Options.log_file)
                fprintf(log_file, '--- new variable solver ends ---\n\n');
                fclose(log_file);
            end

            output = struct;

            % stop the timer
            total_time = toc(total_timer);
            output.total_time = total_time;
        end
    end

    methods(Access = protected)

        function processed_vars = processNewVariables(obj, new_vars)
            % Process the generated new variables to fix potential numerical errors
            % Input:
            %   new_vars: matrix containing the newly generated variables as rows
            % Output:
            %   processed_vars: matrix containing the processed variables as rows

            processed_vars = min(max(new_vars, obj.problem.st2d_constr_lb'), obj.problem.st2d_constr_ub');
        end
    end

    methods(Abstract, Access = protected)

        % Generate the Gurobi MILP model stored in obj.Runtime.Model
        initializeModel(obj);
    end
end

