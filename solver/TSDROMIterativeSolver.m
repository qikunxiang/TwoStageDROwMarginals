classdef TSDROMIterativeSolver < handle
    % Iterative solver that combines the dual LSIP solver and the new variable solver to compute feasible and approximately optimal
    % solutions of the primal and dual formulations of the two-stage distributionally robust optimization problem with marginal
    % constraints
    
    properties(SetAccess = protected, GetAccess = public)
        problem = [];

        DualLSIPSolver;

        NewVariableSolver;

        Options;
    end
    
    methods(Access = public)
        function obj = TSDROMIterativeSolver(DualLSIPSolver, NewVariableSolver, opt)
            % Constructor
            % Inputs:
            %   DualLSIPSolver: object of TSDROMDualLSIPSolver that has not been initialized (the problem property is unset)
            %   NewVariableSolver: object of TSDROMNewVariableSolver that has not been initialized (the problem property is unset)
            %   opt: struct containing the options for the iterative solver
            
            obj.DualLSIPSolver = DualLSIPSolver;
            obj.NewVariableSolver = NewVariableSolver;

            if ~exist('opt', 'var') || isempty(opt)
                opt = struct;
            end

            obj.Options = opt;
    
            if ~isfield(obj.Options, 'display') || isempty(obj.Options.display)
                obj.Options.display = false;
            end

            if ~isfield(obj.Options, 'log_file') || isempty(obj.Options.log_file)
                obj.Options.log_file = '';
            end

            % Whether to recycle active OT constraints
            if ~isfield(obj.Options, 'recycle_OT_constraints') || isempty(obj.Options.recycle_OT_constraints)
                obj.Options.recycle_OT_constraints = true;
            end

            % The algorithm is divided into four phases
            %   I. the exploring phase
            %       This phase is characterized by the dual LSIP solver terminating after a single iteration. In this phase, the goal
            %       is to generate as many dual points as quickly as possible to yield non-trivial dual LSIP solutions. 
            %  II. the squeezing phase
            %       In this phase, the redundant dual points are gradually removed to speed up the dual LSIP solver. At the same time,
            %       the tolerance for upper approximations of the Kantorovich potential functions and the MIP gap are gradually 
            %       decreased to tighten the bounds.
            % III. the stablizing phase that focuses on the lower bound
            %       In this phase, the tolerance for upper approximations of the Kantorovich potential functions and the MIP gap are 
            %       fixed and redundant dual points are no longer removed. An optional time limit can be placed on the MIP solver in
            %       order to generate good candidates of dual points while avoiding wasting time on improving the upper bound.
            %  IV. the stablizing phase that focuses on the upper bound (optional)
            %       In this phase, the tolerance for upper approximations of the Kantorovich potential functions and the MIP gap are 
            %       fixed and redundant dual points are no longer removed. There is also no time limit on the MIP solver in order to
            %       compute a satisfactory upper bound within a handful of iterations.

            % options for phase I
            if ~isfield(obj.Options, 'phase1') || isempty(obj.Options.phase1)
                obj.Options.phase1 = struct;
            end

            % tolerance for upper approximations of the Kantorovich potential functions; specified as a multiplier of the overall
            % tolerance minus the OT tolerance then divided by the number of uncertain variables
            if ~isfield(obj.Options.phase1, 'Kpot_tolerance') || isempty(obj.Options.phase1.Kpot_tolerance)
                obj.Options.phase1.Kpot_tolerance = 100;
            end

            % relative tolerance for MIP gap when solving the global maximization problem
            if ~isfield(obj.Options.phase1, 'MIPGap_tolerance') || isempty(obj.Options.phase1.MIPGap_tolerance)
                obj.Options.phase1.MIPGap_tolerance = 1;
            end

            % threshold for the number of iterations by the dual LSIP solver to switch to phase II
            if ~isfield(obj.Options.phase1, 'DualLSIP_iter_threshold') || isempty(obj.Options.phase1.DualLSIP_iter_threshold)
                obj.Options.phase1.DualLSIP_iter_threshold = 1;
            end

            % options for phase II
            if ~isfield(obj.Options, 'phase2') || isempty(obj.Options.phase2)
                obj.Options.phase2 = struct;
            end

            if ~isfield(obj.Options.phase2, 'DualLSIPOptions') || isempty(obj.Options.phase2.DualLSIPOptions)
                obj.Options.phase2.DualLSIPOptions = struct;
            end

            if ~isfield(obj.Options.phase2, 'DualLSIPLPOptions') || isempty(obj.Options.phase2.DualLSIPLPOptions)
                obj.Options.phase2.DualLSIPLPOptions = struct;
            end

            if ~isfield(obj.Options.phase2, 'DualLSIPGlobalOptions') || isempty(obj.Options.phase2.DualLSIPGlobalOptions)
                obj.Options.phase2.DualLSIPGlobalOptions = struct;
            end

            if ~isfield(obj.Options.phase2, 'NewVarOptions') || isempty(obj.Options.phase2.NewVarOptions)
                obj.Options.phase2.NewVarOptions = struct;
            end

            if ~isfield(obj.Options.phase2, 'NewVarGlobalOptions') || isempty(obj.Options.phase2.NewVarGlobalOptions)
                obj.Options.phase2.NewVarGlobalOptions = struct;
            end

            % threshold for the dual points to be dropped as redundant if not active after a number of iterations
            if ~isfield(obj.Options.phase2, 'dual_points_redundant_threshold') ...
                    || isempty(obj.Options.phase2.dual_points_redundant_threshold)
                obj.Options.phase2.dual_points_redundant_threshold = 10;
            end

            % total number of iterations in phase II
            if ~isfield(obj.Options.phase2, 'iter_num') || isempty(obj.Options.phase2.iter_num)
                obj.Options.phase2.iter_num = 200;
            end

            % tolerance for upper approximations of the Kantorovich potential functions at the beginning of phase II
            if ~isfield(obj.Options.phase2, 'Kpot_tolerance_init') || isempty(obj.Options.phase2.Kpot_tolerance_init)
                obj.Options.phase2.Kpot_tolerance_init = 100;
            end

            % tolerance for upper approximations of the Kantorovich potential functions at the end of phase II
            if ~isfield(obj.Options.phase2, 'Kpot_tolerance_final') || isempty(obj.Options.phase2.Kpot_tolerance_final)
                obj.Options.phase2.Kpot_tolerance_final = 0.5;
            end

            % relative tolerance for MIP gap when solving the global maximization problem at the beginning of phase II
            if ~isfield(obj.Options.phase2, 'MIPGap_tolerance_init') || isempty(obj.Options.phase2.MIPGap_tolerance_init)
                obj.Options.phase2.MIPGap_tolerance_init = 1;
            end

            % relative tolerance for MIP gap when solving the global maximization problem at the end of phase II
            if ~isfield(obj.Options.phase2, 'MIPGap_tolerance_final') || isempty(obj.Options.phase2.MIPGap_tolerance_final)
                obj.Options.phase2.MIPGap_tolerance_final = 0.4;
            end

            % options for phase III
            if ~isfield(obj.Options, 'phase3') || isempty(obj.Options.phase3)
                obj.Options.phase3 = struct;
            end

            if ~isfield(obj.Options.phase3, 'DualLSIPOptions') || isempty(obj.Options.phase3.DualLSIPOptions)
                obj.Options.phase3.DualLSIPOptions = struct;
            end

            if ~isfield(obj.Options.phase3, 'DualLSIPLPOptions') || isempty(obj.Options.phase3.DualLSIPLPOptions)
                obj.Options.phase3.DualLSIPLPOptions = struct;
            end

            if ~isfield(obj.Options.phase3, 'DualLSIPGlobalOptions') || isempty(obj.Options.phase3.DualLSIPGlobalOptions)
                obj.Options.phase3.DualLSIPGlobalOptions = struct;
            end

            if ~isfield(obj.Options.phase3, 'NewVarOptions') || isempty(obj.Options.phase3.NewVarOptions)
                obj.Options.phase3.NewVarOptions = struct;
            end

            if ~isfield(obj.Options.phase3, 'NewVarGlobalOptions') || isempty(obj.Options.phase3.NewVarGlobalOptions)
                obj.Options.phase3.NewVarGlobalOptions = struct;
            end

            % tolerance for upper approximations of the Kantorovich potential functions; specified as a multiplier of the overall
            % tolerance minus the OT tolerance then divided by the number of uncertain variables
            if ~isfield(obj.Options.phase3, 'Kpot_tolerance') || isempty(obj.Options.phase3.Kpot_tolerance)
                obj.Options.phase3.Kpot_tolerance = 0.5;
            end

            % relative tolerance for MIP gap when solving the global maximization problem
            if ~isfield(obj.Options.phase3, 'MIPGap_tolerance') || isempty(obj.Options.phase3.MIPGap_tolerance)
                obj.Options.phase3.MIPGap_tolerance = 0.4;
            end

            % total number of iterations in phase III
            if ~isfield(obj.Options.phase3, 'iter_num') || isempty(obj.Options.phase3.iter_num)
                obj.Options.phase3.iter_num = inf;
            end

            % options for phase IV
            if ~isfield(obj.Options, 'phase4') || isempty(obj.Options.phase4)
                obj.Options.phase4 = struct;
            end

            if ~isfield(obj.Options.phase4, 'DualLSIPOptions') || isempty(obj.Options.phase4.DualLSIPOptions)
                obj.Options.phase4.DualLSIPOptions = struct;
            end

            if ~isfield(obj.Options.phase4, 'DualLSIPLPOptions') || isempty(obj.Options.phase4.DualLSIPLPOptions)
                obj.Options.phase4.DualLSIPLPOptions = struct;
            end

            if ~isfield(obj.Options.phase4, 'DualLSIPGlobalOptions') || isempty(obj.Options.phase4.DualLSIPGlobalOptions)
                obj.Options.phase4.DualLSIPGlobalOptions = struct;
            end

            if ~isfield(obj.Options.phase4, 'NewVarOptions') || isempty(obj.Options.phase4.NewVarOptions)
                obj.Options.phase4.NewVarOptions = struct;
            end

            if ~isfield(obj.Options.phase4, 'NewVarGlobalOptions') || isempty(obj.Options.phase4.NewVarGlobalOptions)
                obj.Options.phase4.NewVarGlobalOptions = struct;
            end

            % tolerance for upper approximations of the Kantorovich potential functions; specified as a multiplier of the overall
            % tolerance minus the OT tolerance then divided by the number of uncertain variables
            if ~isfield(obj.Options.phase4, 'Kpot_tolerance') || isempty(obj.Options.phase4.Kpot_tolerance)
                obj.Options.phase4.Kpot_tolerance = 0.5;
            end

            % relative tolerance for MIP gap when solving the global maximization problem
            if ~isfield(obj.Options.phase4, 'MIPGap_tolerance') || isempty(obj.Options.phase4.MIPGap_tolerance)
                obj.Options.phase4.MIPGap_tolerance = 0.4;
            end
        end
        
        function setProblem(obj, problem)
            % Set the problem to be solved by the solver
            % Input:
            %   problem: object of TSDROMInstance
            
            if ~isempty(obj.problem)
                error('problem has already been set');
            end

            obj.problem = problem;
            obj.DualLSIPSolver.setProblem(problem);
            obj.NewVariableSolver.setProblem(problem);
        end

        function [st1_deci, wc_sampler, primal_sol, dual_sol, DRO_UB, DRO_LB, error_sub, error_prob, output] = ...
                run(obj, tolerance, init_dual_points, OT_tolerance)
            % Run the iterative solver
            % Inputs:
            %   tolerance: the sub-optimality tolerance
            %   init_dual_points: matrix containing the initial dual points as rows (if it is omitted, a single dual point will be
            %   generated)
            %   OT_tolerance: the sub-optimality tolerance when solving the dual LSIP problems (default is equal to 
            %   tolerance / 2 / obj.problem.unc_num)
            % Outputs:
            %   st1_deci: computed first-stage decision
            %   wc_sampler: sampler of the worst-case probability measure associated with the computed first-stage decision
            %   primal_sol: approximately optimal solution of the primal DRO formulation
            %   dual_sol: approximately optimal solution of the dual DRO formulation
            %   DRO_UB: computed upper bound for the objective value of the two-stage DRO problem
            %   DRO_LB: computed lower bound for the objective value of the two-stage DRO problem
            %   error_sub: computed sub-optimality estimate of the primal and dual solutions
            %   error_prob: computed sub-optimality estimate of the worst-case probability measure
            %   output: struct containing additional output information such as timer information

            % start the timer
            total_timer = tic;

            LSIP_total_time = 0;
            LSIP_LP_time = 0;
            LSIP_global_time = 0;
            NewVar_total_time = 0;

            if isempty(obj.problem)
                error('must first set the problem');
            end

            % open the log file
            if ~isempty(obj.Options.log_file)
                log_file = fopen(obj.Options.log_file, 'a');

                if log_file < 0
                    error('cannot open log file');
                end

                fprintf(log_file, '--- iterative solver starts ---\n');
            end

            if ~exist('init_dual_points', 'var') || isempty(init_dual_points)
                init_dual_point = obj.problem.uncertain.computeOneDualPoint();
                dual_points = init_dual_point';
            else
                dual_points = init_dual_points;
            end

            if ~exist('OT_tolerance', 'var') || isempty(OT_tolerance)
                OT_tolerance = tolerance / 2 / obj.problem.unc_num;
            end

            phase = 1;
            tolerance_remainder = tolerance - OT_tolerance * obj.problem.unc_num;
            Kpot_tolerance = obj.Options.phase1.Kpot_tolerance * tolerance_remainder / obj.problem.unc_num;
            MIPGap_tolerance = obj.Options.phase1.MIPGap_tolerance;
            iter = 0;
            DRO_LB = -inf;
            DRO_UB = inf;
            DRO_LB_cand = -inf; %#ok<NASGU>
            DRO_UB_cand = inf; %#ok<NASGU>
            DRO_best_suboptim = inf;
            primal_sol_cand = []; %#ok<NASGU>
            dual_sol_cand = []; %#ok<NASGU>

            stats_size = 100;
            stats = struct( ...
                'support_size', num2cell(zeros(stats_size, 1)), ...
                'DRO_LB', num2cell(zeros(stats_size, 1)), ...
                'DRO_UB', num2cell(zeros(stats_size, 1)), ...
                'error_sub', num2cell(zeros(stats_size, 1)));

            phase2_Kpot_tolerance_scheme = exp(linspace(log(obj.Options.phase2.Kpot_tolerance_init), ...
                log(obj.Options.phase2.Kpot_tolerance_final), obj.Options.phase2.iter_num)') ...
                * tolerance_remainder / obj.problem.unc_num;
            phase2_Kpot_tolerance_scheme = [phase2_Kpot_tolerance_scheme; phase2_Kpot_tolerance_scheme(end)];

            phase2_MIPGap_tolerance_scheme = exp(linspace(log(obj.Options.phase2.MIPGap_tolerance_init), ...
                log(obj.Options.phase2.MIPGap_tolerance_final), obj.Options.phase2.iter_num)');
            phase2_MIPGap_tolerance_scheme = [phase2_MIPGap_tolerance_scheme; phase2_MIPGap_tolerance_scheme(end)];


            while true
                elapsed_time = toc(total_timer);

                % display
                if obj.Options.display
                    fprintf(['Iteration %5d, %5.0f sec, phase %d: ' ...
                        '%6d dual points, LB = %11.4f, UB = %11.4f, sub-optimality = %12.6f, best sub-optimality = %12.6f\n'], ...
                        iter, elapsed_time, phase, size(dual_points, 1), DRO_LB, DRO_UB, DRO_UB - DRO_LB, DRO_best_suboptim);
                end
    
                % logging
                if ~isempty(obj.Options.log_file)
                    fprintf(log_file, ['Iteration %5d, %5.0f sec, phase %d: ' ...
                        '%6d dual points, LB = %11.4f, UB = %11.4f, sub-optimality = %12.6f, best sub-optimality = %12.6f\n'], ...
                        iter, elapsed_time, phase, size(dual_points, 1), DRO_LB, DRO_UB, DRO_UB - DRO_LB, DRO_best_suboptim);
                end

                if DRO_UB - DRO_LB <= tolerance
                    break;
                end

                obj.DualLSIPSolver.setDualPoints(dual_points);

                if phase == 1 || ~obj.Options.recycle_OT_constraints
                    init_constr = obj.DualLSIPSolver.generateInitialOTConstraints();
                else
                    init_constr = obj.DualLSIPSolver.generateInitialOTConstraints([active_dual_probs; ...
                        ones(1, size(dual_points, 1)) / size(dual_points, 1)]);
                end

                LSIP_output = obj.DualLSIPSolver.run(init_constr, OT_tolerance);
                LSIP_total_time = LSIP_total_time + LSIP_output.total_time;
                LSIP_LP_time = LSIP_LP_time + LSIP_output.LP_time;
                LSIP_global_time = LSIP_global_time + LSIP_output.global_time;

                DRO_LB = obj.DualLSIPSolver.getObjectiveLowerBound();

                LSIP_primal_sol = obj.DualLSIPSolver.Runtime.PrimalSolution;
                LSIP_dual_sol = obj.DualLSIPSolver.Runtime.DualSolution;
                LSIP_UB = obj.computePrimalObjective(LSIP_dual_sol);
                active_dual_probs = obj.DualLSIPSolver.getActiveDualPointProbabilities();

                dual_sol = obj.processDualSolution(dual_points, LSIP_primal_sol);

                primal_sol = struct;
                primal_sol.st1_deci = LSIP_dual_sol.st1_deci;
                primal_sol.Kpot_duals_cell = cell(obj.problem.unc_num, 1);

                for marg_id = 1:obj.problem.unc_num
                    primal_sol.Kpot_duals_cell{marg_id} = obj.processKantorovichPotential(LSIP_dual_sol.Kpot_duals_cell{marg_id}, ...
                        Kpot_tolerance);
                end

                obj.NewVariableSolver.setPrimalSolution(primal_sol);
                NewVar_GlobalOptions = struct('MIPGap', MIPGap_tolerance);
                obj.NewVariableSolver.updateOptions(struct, NewVar_GlobalOptions);
                [pool_st2d, vio, NewVar_output] = obj.NewVariableSolver.run();
                NewVar_total_time = NewVar_total_time + NewVar_output.total_time;

                DRO_UB = obj.computePrimalObjective(primal_sol) + vio;

                % decrease each Kantorovich potential function to make it primally feasible
                primal_sol.Kpot_duals_cell{1}(:, 2) = primal_sol.Kpot_duals_cell{1}(:, 2) - vio;

                % store the best candidates of primal and dual solutions encountered so far; these variables are not used in the
                % outputs but are kept in case the algorithm has to be manually terminated prematurely
                if DRO_UB - DRO_LB < DRO_best_suboptim
                    DRO_best_suboptim = DRO_UB - DRO_LB;
                    DRO_UB_cand = DRO_UB; %#ok<NASGU>
                    DRO_LB_cand = DRO_LB;  %#ok<NASGU>
                    primal_sol_cand = primal_sol;  %#ok<NASGU>
                    dual_sol_cand = dual_sol;  %#ok<NASGU>
                end

                iter = iter + 1;

                % statistics of this iteration
                if iter > size(stats, 1)
                    stats = [stats; struct( ...
                        'support_size', num2cell(zeros(stats_size, 1)), ...
                        'DRO_LB', num2cell(zeros(stats_size, 1)), ...
                        'DRO_UB', num2cell(zeros(stats_size, 1)), ...
                        'error_sub', num2cell(zeros(stats_size, 1)))]; %#ok<AGROW>
                    stats_size = stats_size * 2;
                end

                stats(iter).support_size = length(dual_sol.dual_points_probs);
                stats(iter).DRO_LB = DRO_LB;
                stats(iter).DRO_UB = DRO_UB;
                stats(iter).error_sub = DRO_UB - DRO_LB;

                % if the number of iteration in the dual LSIP solver exceed the specified threshold, switch to phase II
                if phase == 1 && LSIP_output.iter > obj.Options.phase1.DualLSIP_iter_threshold
                    phase = 2;

                    obj.DualLSIPSolver.updateOptions(obj.Options.phase2.DualLSIPOptions, ...
                        obj.Options.phase2.DualLSIPLPOptions, obj.Options.phase2.DualLSIPGlobalOptions);
                    obj.NewVariableSolver.updateOptions(obj.Options.phase2.NewVarOptions, obj.Options.phase2.NewVarGlobalOptions);

                    phase2_iter = 0;

                    % start to remove redundant dual points
                    dual_points_age = inf(size(dual_points, 1), 1);
                end

                if phase == 2
                    phase2_iter = phase2_iter + 1;

                    Kpot_tolerance = phase2_Kpot_tolerance_scheme(phase2_iter);
                    MIPGap_tolerance = phase2_MIPGap_tolerance_scheme(phase2_iter);

                    % update the age of the dual points and filter out those too old
                    dual_points_age = dual_points_age + 1;
                    dual_points_age(LSIP_primal_sol.probs_raw > 0) = 0;
                    dual_points_keep_list = dual_points_age <= obj.Options.phase2.dual_points_redundant_threshold;
                    dual_points = [dual_points(dual_points_keep_list, :); pool_st2d];
                    dual_points_age = [dual_points_age(dual_points_keep_list); zeros(size(pool_st2d, 1), 1)];
                    active_dual_probs = [active_dual_probs(:, dual_points_keep_list), ...
                        zeros(size(active_dual_probs, 1), size(pool_st2d, 1))];
                else
                    dual_points = [dual_points; pool_st2d]; %#ok<AGROW>
                    active_dual_probs = [active_dual_probs, zeros(size(active_dual_probs, 1), size(pool_st2d, 1))]; %#ok<AGROW>
                end

                % switch to phase III
                if phase == 2 && phase2_iter > obj.Options.phase2.iter_num
                    phase = 3;

                    phase3_iter = 0;

                    obj.DualLSIPSolver.updateOptions(obj.Options.phase3.DualLSIPOptions, ...
                        obj.Options.phase3.DualLSIPLPOptions, obj.Options.phase3.DualLSIPGlobalOptions);
                    obj.NewVariableSolver.updateOptions(obj.Options.phase3.NewVarOptions, obj.Options.phase3.NewVarGlobalOptions);

                    Kpot_tolerance = obj.Options.phase3.Kpot_tolerance * tolerance_remainder / obj.problem.unc_num;
                    MIPGap_tolerance = obj.Options.phase3.MIPGap_tolerance;
                end

                if phase == 3
                    phase3_iter = phase3_iter + 1;
                end

                % switch to phase IV
                if phase == 3 && phase3_iter > obj.Options.phase3.iter_num
                    phase = 4;

                    obj.DualLSIPSolver.updateOptions(obj.Options.phase4.DualLSIPOptions, ...
                        obj.Options.phase4.DualLSIPLPOptions, obj.Options.phase4.DualLSIPGlobalOptions);
                    obj.NewVariableSolver.updateOptions(obj.Options.phase4.NewVarOptions, obj.Options.phase4.NewVarGlobalOptions);

                    Kpot_tolerance = obj.Options.phase4.Kpot_tolerance * tolerance_remainder / obj.problem.unc_num;
                    MIPGap_tolerance = obj.Options.phase4.MIPGap_tolerance;
                end
            end

            % close the log file
            if ~isempty(obj.Options.log_file)
                fprintf(log_file, '--- iterative solver ends ---\n\n');
                fclose(log_file);
            end

            st1_deci = primal_sol.st1_deci;
            wc_sampler = TSDROMWorstCaseSampler(obj.problem, primal_sol, dual_sol);

            error_sub = DRO_UB - DRO_LB;
            error_prob = LSIP_UB - DRO_LB + vio;

            output = struct;
            output.iter = iter;
            output.LSIP_total_time = LSIP_total_time;
            output.LSIP_LP_time = LSIP_LP_time;
            output.LSIP_global_time = LSIP_global_time;
            output.NewVar_total_time = NewVar_total_time;
            output.stats = stats(1:iter);
            
            % stop the timer
            total_time = toc(total_timer);
            output.total_time = total_time;
        end
    end

    methods(Access = protected)
        function Kpot_updated = processKantorovichPotential(~, Kpot, Kpot_tolerance)
            % Process a convex and piece-wise affine Kantorovich potential function by removing some knots without changing the
            % function values too much
            % Input:
            %   Kpot: matrix with four columns where column 1 represents locations of the knots, column 2 represents the function
            %   values, and columns 3 & 4 represent the left & right derivatives
            %   Kpot_tolerance: tolerance value for the approximation
            % Output:
            %   Kpot_updated: matrix with four columns representing the processed Kantorovich potential function

            knot_num = size(Kpot, 1);

            % with less than three knots there is no need to process the function
            if knot_num <= 2
                Kpot_updated = Kpot;
                return;
            end

            % use two knot indices: the first index corresponds to the current rightmost knot that is kept, the second index
            % corresponds to the current rightmost knot that can be potentially removed

            % the leftmost and rightmost knots are always kept
            keep_index = 1;

            % counter variable that keeps track of the number of knots kept
            knot_counter = 1;

            Kpot_updated = zeros(knot_num, 4);
            Kpot_updated(1, :) = Kpot(1, :);

            while keep_index < knot_num - 1
                % attemp to remove subsequent knots to the right of the current knot until the error is too large
                next_index = keep_index + 1;

                curr_slope = Kpot(keep_index, 4);
                
                while next_index < knot_num
                    % use the knot to the right of next_index to interpolate all knots between keep_index and next_index + 1
                    interp_slope = (Kpot(next_index + 1, 2) - Kpot(keep_index, 2)) / (Kpot(next_index + 1, 1) - Kpot(keep_index, 1));
                    interp_vals =  interp_slope * (Kpot(keep_index + 1:next_index, 1) - Kpot(keep_index, 1)) + Kpot(keep_index, 2);
                    
                    interp_errs = interp_vals - Kpot(keep_index + 1:next_index, 2);

                    if max(interp_errs) <= Kpot_tolerance
                        % maximum error is tolerable, attempt to remove one more knot
                        next_index = next_index + 1;
                        curr_slope = interp_slope;
                    else
                        break;
                    end
                end

                % update the right derivative of the current knot
                Kpot_updated(knot_counter, 4) = curr_slope;

                knot_counter = knot_counter + 1;
                Kpot_updated(knot_counter, 1) = Kpot(next_index, 1);
                Kpot_updated(knot_counter, 2) = Kpot(next_index, 2);
                Kpot_updated(knot_counter, 3) = curr_slope;
                Kpot_updated(knot_counter, 4) = Kpot(next_index, 4);

                keep_index = next_index;
            end

            % always keep the last knot
            if keep_index < knot_num
                knot_counter = knot_counter + 1;
                Kpot_updated(knot_counter, :) = Kpot(end, :);
            end

            Kpot_updated = Kpot_updated(1:knot_counter, :);
        end

        function dual_sol = processDualSolution(~, dual_points, LSIP_sol)
            % Process the dual solution computed by the dual LSIP solver
            % Inputs:
            %   dual_points: matrix containing the dual points as rows
            %   LSIP_sol: struct returned by the dual LSIP solver
            % Output:
            %   dual_sol: struct containing the processed dual solution

            dual_sol = struct;
            dual_sol.dual_in = LSIP_sol.dual_in;
            dual_sol.dual_eq = LSIP_sol.dual_eq;

            pos_list = LSIP_sol.probs > 0;
            dual_sol.dual_points = dual_points(pos_list, :);
            dual_sol.dual_points_probs = LSIP_sol.probs(pos_list);
        end

        function val = computePrimalObjective(obj, primal_sol)
            % Compute the objective value of a (possibly infeasible) primal solution
            % Inputs:
            %   primal_sol: struct representing the primal solution
            % Output:
            %   val: computed objective value

            Kpot_integrals = zeros(obj.problem.unc_num, 1);

            for marg_id = 1:obj.problem.unc_num
                Kpot = primal_sol.Kpot_duals_cell{marg_id};
                marginal = obj.problem.uncertain.marginals{marg_id};

                int_lb = Kpot(:, 3);
                int_ub = [Kpot(1:end - 1, 4); Kpot(end, 4) + 1e-14];
                Kpot_integrals(marg_id) = sum(marginal.computeAffIntegral(Kpot(:, 1), -Kpot(:, 2), int_lb, int_ub));
            end

            val = obj.problem.st1_obj' * primal_sol.st1_deci + sum(Kpot_integrals);
        end
    end
end

