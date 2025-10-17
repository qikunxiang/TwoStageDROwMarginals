classdef (Abstract) LSIPMinCuttingPlaneAlgo < handle
    % Class for cutting-plane algorithms used for solving linear semi-infinite programming (LSIP) problems (minimization)

    properties(GetAccess = public, SetAccess = protected)
        % struct for storing the options for the cutting-plane algorithm
        Options;

        % struct for storing the additional options for the gurobi linear programming (LP) solver
        LPOptions;

        % struct for storing the options for the global optimization oracle
        GlobalOptions;

        % struct for storing pre-computed quantities
        Storage;

        % struct for storing temporary information during runtime
        Runtime;
    end

    methods(Access = public)
        function obj = LSIPMinCuttingPlaneAlgo(opt, LP_opt, gl_opt)
            % Constructor method
            % Inputs:
            %   opt: struct containing the options for the cutting-plane algorithm with the following fields
            %       opt.time_limit: maximum time permitted for the algorithm (default is inf)
            %       opt.display: boolean indicating whether to display information in each iteration (default is true)
            %       opt.log_file: path to the log file where the outputs will be written (default is '')
            %       opt.rescue: struct containing parameters for saving and loading save files
            %       opt.rescue.save_file: path to save the progress (default is '')
            %       opt.rescue.save_interval: interval to save progress (default is 100)
            %       opt.rescue.load_file: path to load a save file before running the algorithm (default is '')
            %   LP_opt: struct for storing the additional options for the gurobi linear programming (LP) solver (default is struct())
            %   gl_opt: struct for storing the options for the global optimization oracle (default is struct())

            % set the default options
            if ~exist('opt', 'var') || isempty(opt)
                opt = struct;
            end

            if ~isfield(opt, 'time_limit') || isempty(opt.time_limit)
                opt.time_limit = inf;
            end

            if ~isfield(opt, 'display') || isempty(opt.display)
                opt.display = true;
            end

            if ~isfield(opt, 'log_file') || isempty(opt.log_file)
                opt.log_file = '';
            end

            if ~isfield(opt, 'rescue') || isempty(opt.rescue)
                opt.rescue = struct;
            end

            if ~isfield(opt.rescue, 'save_file') || isempty(opt.rescue.save_file)
                opt.rescue.save_file = '';
            end

            if ~isfield(opt.rescue, 'save_interval') || isempty(opt.rescue.save_interval)
                opt.rescue.save_interval = 100;
            end

            obj.Options = opt;

            if ~exist('LP_opt', 'var') || isempty(LP_opt)
                LP_opt = struct;
            end

            obj.LPOptions = LP_opt;

            if ~exist('gl_opt', 'var') || isempty(gl_opt)
                gl_opt = struct;
            end

            obj.GlobalOptions = gl_opt;
        end

        function output = run(obj, init_constr, tolerance, load_from_file)
            % Run the cutting-plane algorithm for linear semi-infinite programming. The computed primal and dual solutions as well as 
            % the upper and lower bounds will be stored in the runtime environment afterwards.
            % Inputs:
            %   init_constr: information about the initial constraints
            %   tolerance: numerical tolerance value (default is 1e-4)
            %   load_from_file: string containing the path to the file to load from; the file must be a save from this function
            % Outputs:
            %   output: struct containing additional outputs
            %       output.total_time: total time spent in the algorithm
            %       output.iter: number of iterations
            %       output.LP_time: time spent solving LP problems
            %       output.global_time: time spent in the global minimization oracle

            % start the timer
            total_timer = tic;

            time_limit_exceeded = false;

            if ~exist('tolerance', 'var') || isempty(tolerance)
                tolerance = 1e-4;
            end

            if ~exist('load_from_file', 'var') || isempty(load_from_file)
                loading_from_file = false;
            else
                loading_from_file = true;
            end

            % initialize the algorithm by preparing some constant quantities
            obj.initializeBeforeRun();

            % parameters of the LP solver
            LP_options = struct;

            % disable output by default
            LP_options.OutputFlag = 0;

            % since the gurobi MATLAB interface does not support callback, set a 30-minute time limit on the LP solver; if the solver
            % does not converge when the time limit is hit, it assumes that numerical issues have occurred and restarts the solution
            % process with adjusted parameters for better numerical focus
            LP_options.TimeLimit = 1800;

            % set the additional parameters for the LP solver
            LP_options_fields = fieldnames(obj.LPOptions);
            LP_options_values = struct2cell(obj.LPOptions);

            for fid = 1:length(LP_options_fields)
                LP_options.(LP_options_fields{fid}) = LP_options_values{fid};
            end

            % save the actual options in the object
            obj.LPOptions = LP_options;

            % prepare the runtime environment
            obj.prepareRuntime();

            % generate the inital linear programming minimization problem in gurobi
            obj.Runtime.CurrentLPModel = obj.generateInitialMinModel();

            % initialize the LSIP lower and upper bounds
            obj.Runtime.LSIP_UB = inf;
            obj.Runtime.LSIP_LB = -inf;
            obj.Runtime.iter = 0;

            % indicate that the cutting-plane algorithm is being executed
            obj.Runtime.Finished = false;

            % add the initial constraints to the LP model stored in the runtime environment
            obj.addConstraints(init_constr);

            % some statistics
            LP_time = 0;
            global_time = 0;
            total_time = 0;

            % open the log file
            if ~isempty(obj.Options.log_file)
                log_file = fopen(obj.Options.log_file, 'a');

                if log_file < 0
                    error('cannot open log file');
                end

                fprintf(log_file, '--- cutting-plane algorithm starts ---\n');
            end

            if loading_from_file
                % load all states from the load file
                loaded_save_file = load(load_from_file);
                obj.Runtime = loaded_save_file.saved_runtime;
            end

            % loop until the gap between lower and upper bounds is below the tolerance 
            while true
                % solve the LP problem
                LP_timer = tic;
                LP_trial_num = 0;

                % start with the provided options
                LP_options_runtime = LP_options;
                
                while true
                    LP_result = gurobi(obj.Runtime.CurrentLPModel, LP_options_runtime);

                    if strcmp(LP_result.status, 'OPTIMAL')
                        break;
                    end

                    % if the provided options fail, increase the trial number and generate updated options from the error handler
                    LP_trial_num = LP_trial_num + 1;

                    LP_options_runtime = obj.handleLPErrors(LP_result, LP_options_runtime, LP_trial_num);
                end

                LP_time = LP_time + toc(LP_timer);

                % update the runtime environment
                obj.updateRuntimeAfterLP(LP_result);

                % update the lower bound
                obj.updateLSIPLB(LP_result);

                % call the global minimization oracle
                global_timer = tic;
                [min_lb, optimizers] = obj.callGlobalMinOracle(LP_result.x);
                global_time = global_time + toc(global_timer);

                % update the upper bound
                obj.updateLSIPUB(min_lb, optimizers);

                display_string = obj.buildMessage();

                % display output
                if obj.Options.display
                    fprintf(display_string);
                end

                % write log
                if ~isempty(obj.Options.log_file)
                    fprintf(log_file, display_string);
                end

                % check the termination criterion here
                if obj.checkTerminationCondition(tolerance, min_lb, optimizers)
                    break;
                end

                % reduce cuts
                obj.reduceConstraints(LP_result);

                % generate new constraints
                obj.addConstraints(optimizers);

                % update the iteration counter
                obj.Runtime.iter = obj.Runtime.iter + 1;

                % overwrite manual_save to true to save states while debugging
                manual_save = false;

                % save states
                if ~isempty(obj.Options.rescue.save_file) && (mod(obj.Runtime.iter, obj.Options.rescue.save_interval) == 0 ...
                        || manual_save)
                    saved_runtime = obj.Runtime;
                    save(obj.Options.rescue.save_file, 'saved_runtime');
                end
                
                % check if the elapsed time has exceeded the time limit
                if toc(total_timer) > obj.Options.time_limit
                    warning('LSIP time limit exceeded');

                    if ~isempty(obj.Options.log_file)
                        fprintf(log_file, '%s: LSIP time limit exceeded\n', class(obj));
                    end

                    time_limit_exceeded = true;
                    break;
                end
            end

            % prepare the approximately optimal primal solution
            primal_sol = obj.buildPrimalSolution(LP_result, min_lb);
            obj.Runtime.PrimalSolution = primal_sol;

            % prepare the approximately optimal dual solution
            dual_sol = obj.buildDualSolution(LP_result);
            obj.Runtime.DualSolution = dual_sol;

            % set this flag to indicate that the cutting-plane algorithm has finished execution
            obj.Runtime.Finished = true;

            % close the log file
            if ~isempty(obj.Options.log_file)
                fprintf(log_file, '--- cutting-plane algorithm ends ---\n\n');
                fclose(log_file);
            end

            obj.cleanUpRuntimeAfterRun();

            % prepare additional outputs
            output = struct;
            output.iter = obj.Runtime.iter;
            output.LP_time = LP_time;
            output.global_time = global_time;
            output.time_limit_exceeded = time_limit_exceeded;

            % stop the timer
            total_time = total_time + toc(total_timer);
            output.total_time = total_time;
        end
    end

    methods(Access = protected)
        function prepareRuntime(obj)
            % Prepare the runtime environment by initializing some variables
            obj.Runtime = struct;

            % initialize the warm-start basis used by gurobi
            obj.Runtime.cbasis = [];
            obj.Runtime.vbasis = [];
        end

        function updateRuntimeAfterLP(obj, result)
            % Update the runtime environment after solving each LP
            % Input:
            %   result: struct produced by gurobi

            % update the warm-start basis used by gurobi
            if isfield(result, 'vbasis') && ~isempty(result.vbasis) && isfield(result, 'cbasis') && ~isempty(result.cbasis)
                obj.Runtime.vbasis = result.vbasis;
                obj.Runtime.cbasis = result.cbasis;
            else
                % if warm-start basis information is unavailable from the result, remove these basis to avoid future solvers using
                % stale basis
                obj.Runtime.vbasis = [];
                obj.Runtime.cbasis = [];
            end
        end

        function cleanUpRuntimeAfterRun(obj) %#ok<MANU>
            % Clean up the runtime environment after finishing the cutting-plane algorithm
        end

        function updateLSIPLB(obj, result)
            % Update the LSIP lower bound in the runtime environment after solving each LP 
            % Inputs:
            %   result: struct produced by gurobi
            
            obj.Runtime.LSIP_LB = result.objval;
        end

        function display_string = buildMessage(obj)
            % Build a string as the message to display after each iteration
            % Output:
            %   display_string: the string to display

            display_string = sprintf('%s: iteration %4d: LB = %10.4f, UB = %10.4f, UB - LB = %10.6f\n', ...
                class(obj), ...
                obj.Runtime.iter, ...
                obj.Runtime.LSIP_LB, ...
                obj.Runtime.LSIP_UB, ...
                obj.Runtime.LSIP_UB - obj.Runtime.LSIP_LB);
        end

        function terminate = checkTerminationCondition(obj, tolerance, min_lb, optimizers) %#ok<INUSD>
            % Check whether to terminate the cutting-plane algorithm
            % Inputs:
            %   tolerance: user-specified tolerance value
            %   min_lb: lower bound computed by the global optimization solver
            %   optimizers: computed approximate global optimizers
            % Output:
            %   terminate: boolean value indicating whether to terminate the cutting-plane algorithm

            terminate = obj.Runtime.LSIP_UB - obj.Runtime.LSIP_LB <= tolerance;
        end

        function LP_options_runtime_new = handleLPErrors(obj, LP_result, LP_options_runtime, LP_trial_num) %#ok<INUSD>
            % Handle numerical errors that occurred while solving LP
            % Inputs: 
            %   LP_result: struct returned by the gurobi function representing the result from solving LP
            %   LP_options_runtime: struct containing the current options for solving LP
            %   LP_trial_num: integer representing the number of trials so far
            % Output:
            %   LP_options_runtime_new: struct containing the updated options for solving LP

            LP_options_runtime_new = LP_options_runtime;

            if LP_trial_num == 1
                % if the LP solver has failed once (reaching the time limit without converging), retry after  setting higher numeric 
                % focus, turning off presolve, and removing the existing bases
                LP_options_runtime_new.TimeLimit = LP_options_runtime.TimeLimit * 2;
                LP_options_runtime_new.NumericFocus = 3;
                LP_options_runtime_new.Quad = 1;
                LP_options_runtime_new.Presolve = 0;

                if isfield(obj.Runtime.CurrentLPModel, 'cbasis')
                    obj.Runtime.CurrentLPModel = rmfield(obj.Runtime.CurrentLPModel, 'cbasis');
                end

                if isfield(obj.Runtime.CurrentLPModel, 'vbasis')
                    obj.Runtime.CurrentLPModel = rmfield(obj.Runtime.CurrentLPModel, 'vbasis');
                end
            else
                error('error encountered while solving LP');
            end
        end
    end

    methods(Abstract, Access = protected)
        % Initialize the algorithm by preparing some constant quantities
        initializeBeforeRun(obj);

        % Generate the initial linear programming model for gurobi
        model = generateInitialMinModel(obj);

        % Given a decision vector, call the global minimization oracle to approximately determine the "most violated" constraints and
        % return a lower bound for the minimal value
        [min_lb, optimizers] = callGlobalMinOracle(obj, vec);

        % Given a collection of approximate optimizers from the global minimization oracle, generate and add the corresponding linear
        % constraints
        addConstraints(obj, optimizers);

        % Remove some of the constraints to speed up the LP solver
        reduceConstraints(obj, result);

        % Update the LSIP upper bound after each call to the global minimization oracle
        updateLSIPUB(obj, min_lb, optimizers);

        % Given the output from gurobi and a lower bound for the optimal value of the global minimization oracle, build the 
        % corresponding primal solution
        primal_sol = buildPrimalSolution(obj, result, violation);

        % Given the output from gurobi, build the corresponding dual solution
        dual_sol = buildDualSolution(obj, result);
    end
end

