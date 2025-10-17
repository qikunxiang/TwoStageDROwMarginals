classdef TSDROMDualLSIPSolver < LSIPMinCuttingPlaneAlgo
    % Cutting-plane algorithm for solving the restricted dual formulation of the two-stage distributionally robust optimization problem
    % with marginal constraints

    properties(Constant)
        % Constant added to every element before determining the approximately unique elements
        UNIQUENESS_SHIFT = exp(1);

        % The number of digits to round each element to determine the approximately unique elements
        UNIQUENESS_ROUNDING = 5;
    end
    
    properties(SetAccess = protected, GetAccess = public)
        % The problem instance of class TSDROMInstance
        problem = [];

        % Number of points in the finite collection of the feasible set of the dualized inner minimization problem
        dual_point_num = 0;

        % Matrix containing a finite subcollection of the feasible set of the dualized inner minimization problem as rows
        dual_points = [];

        % Cell array containing the projections of the dual points via the columns of problem.st2d_obj_unc
        proj_dual_points = {};
        
        % Cell array containing the maximum distances between the projected dual points and the supports of the marginals
        proj_max_distances = {};

        % Cell array where each cell is a vector containing the indices mapping the original dual points to the projections
        proj_mappings = {};

        % Cell array where each cell is a vector containing the indices of the projections in the sequence of original dual points from
        % left to right
        proj_indices = {};
    end
    
    methods
        function obj = TSDROMDualLSIPSolver(opt, LP_opt, gl_opt)
            % Constructor
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
            
            obj@LSIPMinCuttingPlaneAlgo(opt, LP_opt, gl_opt);

            % set the default options for reducing constraints
            if ~isfield(obj.Options, 'reduce') || isempty(obj.Options.reduce)
                obj.Options.reduce = struct;
            end

            if ~isfield(obj.Options.reduce, 'thres') || isempty(obj.Options.reduce.thres)
                obj.Options.reduce.thres = inf;
            end

            if ~isfield(obj.Options.reduce, 'freq') || isempty(obj.Options.reduce.freq)
                obj.Options.reduce.freq = 20;
            end

            if ~isfield(obj.Options.reduce, 'max_iter') || isempty(obj.Options.reduce.max_iter)
                obj.Options.reduce.max_iter = inf;
            end

            % set the default options for reducing constraints; two thresholds are used for determining which constraints to remove 
            % based on their slackness values: whenever the (negative) slackness is above obj.Options.reduce.thres or its quantile
            % among all the non-tight constraints is below obj.Options.reduce.thres_quantile, the constraint is removed
            if ~isfield(obj.Options.reduce, 'thres_quantile') || isempty(obj.Options.reduce.thres_quantile)
                obj.Options.reduce.thres_quantile = 1;
            end

            % minimum value of slackness (in absolute value) for a constraint to be considered not tight where a constraint is only
            % flagged as removable only if |slack| > obj.Options.reduce.min_slack
            if ~isfield(obj.Options.reduce, 'min_slack') || isempty(obj.Options.reduce.min_slack)
                obj.Options.reduce.min_slack = 0;
            end

            % boolean indicating whether the initial constraints should be preserved throughout the cutting-plane algorithm
            if ~isfield(obj.Options.reduce, 'preserve_init_constr') || isempty(obj.Options.reduce.preserve_init_constr)
                obj.Options.reduce.preserve_init_constr = true;
            end

            % set the default option for the feasibility tolerance of dual points
            if ~isfield(obj.Options, 'dual_feasibility_tolerance') || isempty(obj.Options.dual_feasibility_tolerance)
                obj.Options.dual_feasibility_tolerance = 1e-8;
            end

            % set the default option for the minimum probability value for an atom in the primal LSIP solution to be kept in the sense
            % that an atom is kept in the primal LSIP solution only if probability > obj.Options.primal_prob_thres
            if ~isfield(obj.Options, 'primal_prob_thres') || isempty(obj.Options.primal_prob_thres)
                obj.Options.primal_prob_thres = 0;
            end

            % set the default option for the minimum contribution value for an atom in a feasible LSIP solution to be kept when solving
            % optimal transport problems for generating new constraints
            if ~isfield(obj.Options, 'OT_contrib_thres') || isempty(obj.Options.OT_contrib_thres)
                obj.Options.OT_contrib_thres = 1e-7;
            end

            % set the default option for the minimum weight value for a Kantorovich potential function in the dual LSIP solution to be 
            % considered in the sense that a function is considered in the dual LSIP solution only if 
            % weight > obj.Options.dual_weight_thres
            if ~isfield(obj.Options, 'dual_weight_thres') || isempty(obj.Options.dual_weight_thres)
                obj.Options.dual_weight_thres = 1e-6;
            end

            % this value controls the threshold below which an approximate optimizer will be used to generate new cuts; using 0 as the
            % threshold might result in a large number of new cuts that are close to being non-violations
            if ~isfield(obj.GlobalOptions, 'objective_threshold') || isempty(obj.GlobalOptions.objective_threshold)
                obj.GlobalOptions.objective_threshold = -5e-7;
            end

            if ~isfield(obj.GlobalOptions, 'display') || isempty(obj.GlobalOptions.display)
                obj.GlobalOptions.display = true;
            end

            if ~isfield(obj.GlobalOptions, 'log_file') || isempty(obj.GlobalOptions.log_file)
                obj.GlobalOptions.log_file = '';
            end
        end

        function updateOptions(obj, opt, LP_opt, gl_opt)
            % Update the options for the solver, the LP solver, and the global solver
            % Inputs:
            %   opt: struct containing the additional options for the solver
            %   LP_opt: struct containing the additional options for the LP solver
            %   gl_opt: struct containing the additional options for the global solver
            
            opt_fields = fieldnames(opt);
            opt_values = struct2cell(opt);

            for fid = 1:length(opt_fields)
                obj.Options.(opt_fields{fid}) = opt_values{fid};
            end

            LP_opt_fields = fieldnames(LP_opt);
            LP_opt_values = struct2cell(LP_opt);

            for fid = 1:length(LP_opt_fields)
                obj.LPOptions.(LP_opt_fields{fid}) = LP_opt_values{fid};
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
            %   problem: object of TSDROMInstance

            if ~isempty(obj.problem)
                error('problem has already been set');
            end

            obj.problem = problem;

            obj.dual_point_num = 0;
            obj.dual_points = [];
            obj.proj_dual_points = {};
            obj.proj_max_distances = {};
            obj.proj_mappings = {};
            obj.proj_indices = {};
        end

        function setDualPoints(obj, points)
            % Set the dual points
            % Input:
            %   points: matrix containing a finite subcollection of the feasible set of the dualized inner maximization problem as rows

            if isempty(obj.problem)
                error('problem must be set first');
            end

            assert(size(points, 2) == obj.problem.st2d_deci_num, 'dimensionality of dual points is incorrect');
            assert(all(all(points' >= obj.problem.st2d_constr_lb)), 'dual points do not satisfy lower bound constraints');
            assert(all(all(points' <= obj.problem.st2d_constr_ub)), 'dual points do not satisfy upper bound constraints');
            assert(all(all(obj.problem.st2d_constr_mat_in * points' - obj.problem.st2d_constr_rhs_in ...
                < obj.Options.dual_feasibility_tolerance)), 'dual points do not satisfy the dual feasibility constraints');

            if ~obj.problem.uncertain.checkDualFeasibility(points)
                error('not enough dual points to make the restricted dual problem feasible');
            end

            obj.dual_point_num = size(points, 1);
            obj.dual_points = points;
            obj.initializeDualPoints();
        end

        function optimizers = generateInitialOTConstraints(obj, probs)
            % Generate a set of initial OT-related inequality constraints
            % Input:
            %   probs: (optional) matrix where each row is a reference probability used to generate pairs of Kantorovich potential 
            %   functions for bounding the OT costs from below (default value is the uniform distribution on all the dual points)
            % Output:
            %   optimizers: same format as the output of the function callGlobalMinOracle

            if isempty(obj.problem)
                error('problem must be set first');
            end

            if isempty(obj.dual_points)
                error('dual points must be set first');
            end

            if ~exist('probs', 'var') || isempty(probs)
                % take the uniform distribution over the dual points
                probs = ones(obj.dual_point_num, 1)';
                probs = probs / sum(probs);
            end

            prob_num = size(probs, 1);

            Kpot_integrals_cell = cell(prob_num, 1);
            Kpot_duals_projeval_cell = cell(prob_num, 1);
            Kpot_duals_cell = cell(prob_num, 1);

            for prob_id = 1:prob_num
                % generate a set of inequality constraints that are sufficient to guarantee that the LP relaxation of the LSIP 
                % minimization problem is bounded from below
                [~, Kpot_integrals_cell{prob_id}, Kpot_duals_projeval_cell{prob_id}, Kpot_duals_cell{prob_id}] = ...
                    obj.computeProjectedOT(probs(prob_id, :)');
            end

            optimizers = struct;
            optimizers.marg_ids = repmat((1:obj.problem.unc_num)', prob_num, 1);
            optimizers.Kpot_integrals = vertcat(Kpot_integrals_cell{:});
            optimizers.Kpot_duals_projeval = vertcat(Kpot_duals_projeval_cell{:});
            optimizers.Kpot_duals_cell = vertcat(Kpot_duals_cell{:});
            optimizers.probs = repelem(probs, obj.problem.unc_num, 1);
        end

        function lb = getObjectiveLowerBound(obj)
            % Get the computed lower bound for the two-stage DRO problem
            % Output:
            %   lb: value of the lower bound

            if ~isfield(obj.Runtime, 'DualSolution') ...
                    || isempty(obj.Runtime.DualSolution)
                error('need to first execute the cutting-plane algorithm');
            end

            lb = -obj.Runtime.LSIP_UB;
        end

        function probs = getActiveDualPointProbabilities(obj)
            % Get a list of probabilities over the dual points which were used at some point of the cutting-plane algorithm to generate
            % OT-related inequality constraints that are active in the computed approximate dual optimizer
            % Output:
            %   probs: matrix where each row represents a probability over the dual points

            if ~obj.Runtime.Finished
                error('This function can only be called once the cutting-plane algorithm has finished execution');
            end

            probs = obj.Runtime.probs(obj.Runtime.active_constr_indices, :);
            [~, ~, umap] = unique(round(probs, 6), 'rows', 'stable');
            probs = probs(umap, :);
        end
    end

    methods(Access = protected)

        function initializeDualPoints(obj)
            % Initialize the dual points by calculating their projections and establish the mappings to the projections

            obj.proj_dual_points = cell(obj.problem.unc_num, 1);
            obj.proj_mappings = cell(obj.problem.unc_num, 1);
            obj.proj_indices = cell(obj.problem.unc_num, 1);

            for marg_id = 1:obj.problem.unc_num
                proj_i = full(obj.dual_points * obj.problem.st2d_obj_unc(:, marg_id));

                [~, indices, mapping] = unique(round(proj_i + obj.UNIQUENESS_SHIFT, obj.UNIQUENESS_ROUNDING), 'sorted');
                obj.proj_dual_points{marg_id} = proj_i(indices);
                obj.proj_max_distances{marg_id} = max(abs(obj.proj_dual_points{marg_id} ...
                    - obj.problem.uncertain.marginal_bounds(marg_id, :)), [], 2);
                obj.proj_mappings{marg_id} = mapping;
                obj.proj_indices{marg_id} = indices;
            end
        end

        function proj_probs = computeProjectedProbs(obj, probs)
            % Given a probability on the finite collection of dual points, compute its projections 
            % Input:
            %   probs: vector containing the probabilities of the dual points
            % Output:
            %   proj_probs: cell array containing the projected probabilities

            proj_probs = cell(obj.problem.unc_num, 1);

            for marg_id = 1:obj.problem.unc_num
                proj_probs{marg_id} = accumarray(obj.proj_mappings{marg_id}, probs, size(obj.proj_dual_points{marg_id}));
                proj_probs{marg_id} = proj_probs{marg_id} / sum(proj_probs{marg_id});
            end
        end

        function [OTcosts, Kpot_integrals, Kpot_duals_projeval, Kpot_duals_cell] = computeProjectedOT(obj, probs)
            % Given a probability on the finite collection of dual points, solve one-dimensional semi-discrete optimal transport of its
            % projections to the marginals
            % Input:
            %   probs: vector containing the probabilities of the dual points
            % Outputs: 
            %   OTcosts: vector containing the optimal transport costs of the one-dimensional OT problems
            %   Kpot_integrals: vector containing the integrals of the Kantorovich potential functions \psi_i with respect to the
            %   marginals
            %   Kpot_duals_projeval: matrix where each row contains the Kantorovich potential function f_i evaluated at the projection
            %   of the dual points
            %   Kpot_duals_cell: cell array where each cell contains the Kantorovich potential function f_i evaluated at the knots from
            %   left to right

            OTcosts = zeros(obj.problem.unc_num, 1);
            Kpot_integrals = zeros(obj.problem.unc_num, 1);
            Kpot_duals_projeval = zeros(obj.problem.unc_num, obj.dual_point_num);
            Kpot_duals_cell = cell(obj.problem.unc_num, 1);
            proj_probs = obj.computeProjectedProbs(probs);

            for marg_id = 1:obj.problem.unc_num
                % filter out atoms that contribute negligibly to the OT cost
                OT_max_contributions = proj_probs{marg_id}' * obj.proj_max_distances{marg_id};
                proj_probs{marg_id}(OT_max_contributions < obj.Options.OT_contrib_thres) = 0;
                proj_probs{marg_id} = proj_probs{marg_id} / sum(proj_probs{marg_id});

                [OTcost, Kpot_integral, Kpot_dual] = ...
                    obj.problem.uncertain.marginals{marg_id}.computeOT(obj.proj_dual_points{marg_id}, proj_probs{marg_id});

                OTcosts(marg_id) = OTcost;
                Kpot_integrals(marg_id) = Kpot_integral;

                % Kpot_dual contains the values of the Kantorovich potential function evaluated at the knots from left to right; we
                % utilize the projection mapping to obtain the values of the Kantorovich potential function evaluated at the
                % projections of the dual points in the order of the dual points
                Kpot_duals_projeval(marg_id, :) = Kpot_dual(obj.proj_mappings{marg_id}, 2);

                Kpot_duals_cell{marg_id} = Kpot_dual;
            end
        end

        function initializeBeforeRun(~)
            % Initialize the algorithm by preparing some constant quantities
            % Does not have to do anything
        end

        function prepareRuntime(obj)
            % Prepare the runtime environment by initializing some variables

            prepareRuntime@LSIPMinCuttingPlaneAlgo(obj);

            obj.Runtime.probs = zeros(0, obj.dual_point_num);
        end
        
        function model = generateInitialMinModel(obj)
            % Generate the initial linear programming model for gurobi
            % Output:
            %   model: struct containing the gurobi model

            % the decision variables:
            % 1. the auxiliary variables encoding the OT costs (length = obj.problem.unc_num)
            % 2. the probabilities of the dual points (length = obj.problem.uncertain.dual_point_num)
            % 3. the dual variables corresponding to the first-stage inequality constraints (length = obj.problem.st1_constr_num_in)
            % 4. the dual variables corresponding to the first-stage equality constraints (length = st1_constr_num_eq)

            model = struct;
            model.modelsense = 'min';
            model.objcon = 0;

            model.indices = struct;
            model.indices.OTcosts = (1:obj.problem.unc_num)';
            counter = obj.problem.unc_num;
            model.indices.probs = counter + (1:obj.dual_point_num)';
            counter = counter + obj.dual_point_num;
            model.indices.dual_in = counter + (1:obj.problem.st1_constr_num_in)';
            counter = counter + obj.problem.st1_constr_num_in;
            model.indices.dual_eq = counter + (1:obj.problem.st1_constr_num_eq)';
            deci_num = counter + obj.problem.st1_constr_num_eq;

            model.obj = [ones(obj.problem.unc_num, 1); ...
                -obj.dual_points * obj.problem.st2d_obj_itc; ...
                -obj.problem.st1_constr_rhs_in; ...
                -obj.problem.st1_constr_rhs_eq];
            model.lb = [-inf(obj.problem.unc_num, 1); ...
                zeros(obj.dual_point_num, 1); ...
                -inf(obj.problem.st1_constr_num_in, 1); ...
                -inf(obj.problem.st1_constr_num_eq, 1)];
            model.ub = [inf(obj.problem.unc_num, 1); ...
                inf(obj.dual_point_num, 1); ...
                zeros(obj.problem.st1_constr_num_in, 1); ...
                inf(obj.problem.st1_constr_num_eq, 1)];
            
            model.constr_indices = struct;

            % inequality constraint about convex combination of the dual points
            A_in = sparse([zeros(obj.problem.st1_deci_num, obj.problem.unc_num), ...
                obj.problem.st2d_obj_act' * obj.dual_points', ...
                -obj.problem.st1_constr_mat_in', ...
                -obj.problem.st1_constr_mat_eq']);
            rhs_in = -obj.problem.st1_obj;
            model.constr_indices.st1 = (1:obj.problem.st1_deci_num)';

            % equality constraint requiring all weights to sum to 1
            A_eq = sparse(ones(obj.dual_point_num, 1), model.indices.probs, 1, 1, deci_num);
            rhs_eq = 1;
            model.constr_indices.sum = obj.problem.st1_deci_num + 1;

            model.A = [A_in; A_eq];
            model.rhs = [rhs_in; rhs_eq];
            model.sense = [repmat('>', obj.problem.st1_deci_num, 1); '='];
            model.init_constr_num = size(model.A, 1);

            % vector for storing the associating between the inequality constraints and the marginals
            model.constr_OT_marg_ids = zeros(0, 1);

            % cell array for storing the values of the Kantorovich potential functions f_i evaluated at the knots from left to right
            model.Kpot_duals_cell = cell(obj.problem.unc_num, 1);

            for marg_id = 1:obj.problem.unc_num
                model.Kpot_duals_cell{marg_id} = zeros(length(obj.proj_dual_points{marg_id}), 4, 0);
            end
        end

        function [min_lb, optimizers] = callGlobalMinOracle(obj, vec)
            % Given a decision vector, call the global minimization oracle to approximately determine the "most violated" constraints 
            % and return a lower bound for the minimal value
            % Input:
            %   vec: vector containing the current LP optimizer
            % Outputs:
            %   min_lb: vector containing the lower bounds for the global minimization problems
            %   optimizers: struct containing information about the cuts to generate

            indices = obj.Runtime.CurrentLPModel.indices;
            OTcosts_lb = vec(indices.OTcosts);

            % the raw probabilities directly extracted from the dual LP optimizer where atoms with negligible probabilities have not
            % been filtered out; the negligible atoms will be filtered out inside obj.computeProjectedOT
            probs_raw = max(vec(indices.probs), 0);

            [OTcosts, Kpot_integrals, Kpot_duals_projeval, Kpot_duals_cell] = obj.computeProjectedOT(probs_raw);

            min_lb = OTcosts_lb - OTcosts;

            optimizers = struct;
            optimizers.marg_ids = (1:obj.problem.unc_num)';
            optimizers.Kpot_integrals = Kpot_integrals;
            optimizers.Kpot_duals_projeval = Kpot_duals_projeval;
            optimizers.Kpot_duals_cell = Kpot_duals_cell;
            optimizers.probs = repelem(probs_raw', obj.problem.unc_num, 1);
        end

        function addConstraints(obj, optimizers)
            % Given a collection of approximate optimizers from the global minimization oracle, generate and add the corresponding 
            % linear constraints
            % Input:
            %   optimizers: struct returned by the function callGlobalMinOracle

            constr_num = length(optimizers.marg_ids);
            indices = obj.Runtime.CurrentLPModel.indices;
            
            A_new = sparse( ...
                repmat((1:constr_num)', 1 + obj.dual_point_num, 1), ...
                [indices.OTcosts(optimizers.marg_ids); ...
                    repelem(indices.probs(1:obj.dual_point_num), constr_num, 1)], ...
                [ones(constr_num, 1); optimizers.Kpot_duals_projeval(:)], ...
                constr_num, length(obj.Runtime.CurrentLPModel.obj));

            rhs_new = -optimizers.Kpot_integrals;
            
            % add the newly generated constraints to the end
            obj.Runtime.CurrentLPModel.A = [obj.Runtime.CurrentLPModel.A; A_new];
            obj.Runtime.CurrentLPModel.rhs = [obj.Runtime.CurrentLPModel.rhs; rhs_new];
            obj.Runtime.CurrentLPModel.sense = [obj.Runtime.CurrentLPModel.sense; repmat('>', constr_num, 1)];
            obj.Runtime.CurrentLPModel.constr_OT_marg_ids = [obj.Runtime.CurrentLPModel.constr_OT_marg_ids; optimizers.marg_ids];
            obj.Runtime.probs = [obj.Runtime.probs; optimizers.probs];

            % store the Kantorovich potentials which will be later used to construct an approximately optimal dual LSIP solution
            for marg_id = 1:obj.problem.unc_num
                obj.Runtime.CurrentLPModel.Kpot_duals_cell{marg_id} = cat(3, obj.Runtime.CurrentLPModel.Kpot_duals_cell{marg_id}, ...
                    optimizers.Kpot_duals_cell{optimizers.marg_ids == marg_id});
            end

            if ~isempty(obj.Runtime.vbasis) && ~isempty(obj.Runtime.cbasis)
                obj.Runtime.CurrentLPModel.vbasis = obj.Runtime.vbasis;
                obj.Runtime.CurrentLPModel.cbasis = [obj.Runtime.cbasis; zeros(constr_num, 1)];
            else
                if isfield(obj.Runtime.CurrentLPModel, 'vbasis')
                    obj.Runtime.CurrentLPModel = rmfield(obj.Runtime.CurrentLPModel, 'vbasis');
                end

                if isfield(obj.Runtime.CurrentLPModel, 'cbasis')
                    obj.Runtime.CurrentLPModel = rmfield(obj.Runtime.CurrentLPModel, 'cbasis');
                end
            end

            % if obj.Runtime.CurrentLPModel.init_OT_constr_num is not set, it means that this is the first call to obj.addConstraints 
            % which generates the initial OT-related constraints; this number is stored in the runtime environment
            if ~isfield(obj.Runtime.CurrentLPModel, 'init_OT_constr_num_in') || isempty(obj.Runtime.CurrentLPModel.init_OT_constr_num)
                obj.Runtime.CurrentLPModel.init_OT_constr_num = constr_num;
            end
        end

        function reduceConstraints(obj, result)
            % Remove some of the constraints to speed up the LP solver
            % Input:
            %   result: output of the gurobi solver

            if isinf(obj.Options.reduce.thres) || obj.Runtime.iter <= 0 || obj.Runtime.iter > obj.Options.reduce.max_iter ...
                    || mod(obj.Runtime.iter, obj.Options.reduce.freq) ~= 0
                return;
            end

            cut_slack_OT = result.slack(obj.Runtime.CurrentLPModel.init_constr_num + 1:end);
            flagged_OT_indices = cut_slack_OT < -obj.Options.reduce.min_slack;

            if obj.Options.reduce.preserve_init_constr && isfield(obj.Runtime.CurrentLPModel, 'init_OT_constr_num') ...
                    && ~isempty(obj.Runtime.CurrentLPModel.init_OT_constr_num)
                flagged_OT_indices(1:obj.Runtime.CurrentLPModel.init_OT_constr_num) = false;
            end

            cut_slack_thres_quantile = quantile(-cut_slack_OT(flagged_OT_indices), obj.Options.reduce.thres_quantile);
            keep_list_OT = ~flagged_OT_indices | (cut_slack_OT >= -obj.Options.reduce.thres ...
                & cut_slack_OT >= -cut_slack_thres_quantile);

            % update the constraints
            keep_list_full = [true(obj.Runtime.CurrentLPModel.init_constr_num, 1); keep_list_OT];
            obj.Runtime.CurrentLPModel.A = obj.Runtime.CurrentLPModel.A(keep_list_full, :);
            obj.Runtime.CurrentLPModel.rhs = obj.Runtime.CurrentLPModel.rhs(keep_list_full);
            obj.Runtime.CurrentLPModel.sense = obj.Runtime.CurrentLPModel.sense(keep_list_full);
            obj.Runtime.probs = obj.Runtime.probs(keep_list_full, :);

            for marg_id = 1:obj.problem.unc_num
                keep_list_marg = keep_list_OT(obj.Runtime.CurrentLPModel.constr_OT_marg_ids == marg_id);
                obj.Runtime.CurrentLPModel.Kpot_duals_cell{marg_id} = ...
                    obj.Runtime.CurrentLPModel.Kpot_duals_cell{marg_id}(:, :, keep_list_marg);
            end

            obj.Runtime.CurrentLPModel.constr_OT_marg_ids = obj.Runtime.CurrentLPModel.constr_OT_marg_ids(keep_list_OT);

            if ~isempty(obj.Runtime.cbasis)
                obj.Runtime.cbasis = obj.Runtime.cbasis(keep_list_full);
            end
        end

        function updateLSIPUB(obj, min_lb, ~)
            % Update the LSIP upper bound after each call to the global minimization oracle
            % Inputs:
            %   min_lb: vector containing the lower bounds for the global minimization problems
            %   optimizers: struct returned by the function callGlobalMinOracle

            obj.Runtime.LSIP_UB = obj.Runtime.LSIP_LB - min(sum(min_lb), 0);
        end

        function primal_sol = buildPrimalSolution(obj, result, min_lb)
            % Given the output from gurobi and a lower bound for the optimal value of the global minimization oracle, build the 
            % corresponding primal solution
            % Inputs:
            %   result: output of the gurobi solver
            %   min_lb: vector containing the lower bounds for the global minimization problems
            % Output:
            %   primal_sol: struct representing the computed primal LSIP solution

            indices = obj.Runtime.CurrentLPModel.indices;

            primal_sol = struct;
            primal_sol.OTcosts = result.x(indices.OTcosts) - min_lb;
            primal_sol.dual_in = result.x(indices.dual_in);
            primal_sol.dual_eq = result.x(indices.dual_eq);

            % the raw probabilities directly extracted from the dual LP optimizer where atoms with negligible probabilities have not
            % been filtered out
            probs_raw = result.x(indices.probs);

            % filter out the probabilities that are negligible, i.e., below obj.Options.primal_prob_thres
            probs_negligible_list = probs_raw < obj.Options.primal_prob_thres;
            probs_keep_indices = find(~probs_negligible_list);
            probs = probs_raw;
            probs_negligible = probs(probs_negligible_list, :);
            dist_mat = pdist2(obj.dual_points(probs_negligible_list, :), obj.dual_points(~probs_negligible_list, :));
            [~, min_indices] = min(dist_mat, [], 2);
            probs_negligible_agg = accumarray(min_indices, probs_negligible, [length(probs_keep_indices), 1]);
            probs(probs_negligible_list) = 0;
            probs(probs_keep_indices) = probs(probs_keep_indices) + probs_negligible_agg;
            probs = probs / sum(probs);
            
            % also compute the projection of the probability measures
            proj_probs = obj.computeProjectedProbs(probs);

            primal_sol.probs_raw = probs_raw;
            primal_sol.probs = probs;
            primal_sol.proj_probs = proj_probs;
        end

        function dual_sol = buildDualSolution(obj, result)
            % Given the output from gurobi, build the corresponding dual solution
            % Input:
            %   result: output of the gurobi solver
            % Output:
            %   dual_sol: struct representing the computed dual LSIP solution
            
            constr_indices = obj.Runtime.CurrentLPModel.constr_indices;

            % the dual variable corresponding to the equality constraint related to the first-stage objective vector is precisely the
            % first-stage decision variable
            st1_deci = result.pi(constr_indices.st1);

            % the dual variable corresponding to the equality constraint requiring the probabilities to sum to 1 can be added to the
            % first Kantorovich potential function f_1
            addconst = result.pi(constr_indices.sum);

            % the dual variables corresponding to the inequality constraints are weights combining the Kantorovich potential functions
            % generated throughout the cutting-plane algorithm;
            % in theory, for each marginal, exactly one Kantorovich potential function f_i will have weight 1 and all the rest of the
            % Kantorovich potential functions will have weight 0;
            % however, to handle the edge case when the weights are split into several Kantorovich potential functions (which might
            % never occur in practice), we compute the convex combination of all the Kantorovich potential functions with these weights
            % for each marginal
            func_combweights = result.pi(obj.Runtime.CurrentLPModel.init_constr_num + 1:end);
            
            % save the indices of the active constraints in the runtime
            obj.Runtime.active_constr_indices = func_combweights >= obj.Options.dual_weight_thres;

            marg_ids = obj.Runtime.CurrentLPModel.constr_OT_marg_ids;
            Kpot_duals_cell = obj.Runtime.CurrentLPModel.Kpot_duals_cell;
            Kpot_duals_opt_cell = cell(obj.problem.unc_num, 1);

            for marg_id = 1:obj.problem.unc_num
                marg_list = marg_ids == marg_id;
                func_combweights_marg = func_combweights(marg_list, :);

                % remove weights that fall below the threshold
                func_combweights_marg(func_combweights_marg < obj.Options.dual_weight_thres) = 0;
                func_combweights_marg = func_combweights_marg / sum(func_combweights_marg);

                Kpot_duals_opt_cell{marg_id} = sum(Kpot_duals_cell{marg_id} .* reshape(func_combweights_marg, 1, 1, []), 3);
            end

            Kpot_duals_opt_cell{1}(:, 2) = Kpot_duals_opt_cell{1}(:, 2) + addconst;

            dual_sol = struct;
            dual_sol.st1_deci = st1_deci;
            dual_sol.Kpot_duals_cell = Kpot_duals_opt_cell;
        end

        function terminate = checkTerminationCondition(obj, OT_tolerance, min_lb, optimizers) %#ok<INUSD>
            % Check whether to terminate the cutting-plane algorithm
            % Inputs:
            %   tolerance: user-specified tolerance value related to discrepancies in OT costs
            %   min_lb: lower bound computed by the global optimization solver
            %   optimizers: computed approximate global optimizers
            % Output:
            %   terminate: boolean value indicating whether to terminate the cutting-plane algorithm

            terminate = -min(min_lb) <= OT_tolerance;
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
            elseif LP_trial_num == 2 && isfield(LP_options_runtime, 'Method') && LP_options_runtime.Method == 2 ...
                && isfield(LP_options_runtime, 'Crossover') && LP_options_runtime.Crossover == 0
                % if the solver is using the barrier algorithm with crossover disabled, switch to the dual-simplex algorithm
                LP_options_runtime_new.Method = 1;
                LP_options_runtime_new = rmfield(LP_options_runtime_new, 'Crossover');
            elseif LP_trial_num == 3 && isfield(LP_options_runtime, 'Method') && LP_options_runtime.Method == 1
                % if the dual-simplex algorithm also fails, switch to the primal-simplex algorithm
                LP_options_runtime_new.Method = 0;
            else
                while true
                    % do nothing
                    warning('waiting for intervention...');
                    pause(10);
                    
                    continue_flag = false;

                    if continue_flag
                        break; %#ok<UNRCH>
                    end
                end
            end
        end
    end
end

