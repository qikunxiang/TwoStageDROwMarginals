function [LSIP_UB, LSIP_LB, coup_meas, ...
    stage1_decision, stage1_multipliers, dual_func, output, options] ...
    = DRO_cutplane(DRO, knots, expval, coup_init, options)
% Solve a two-stage distributionally robust optimization (DRO) problem with
% marginal constraints via the cutting plane algorithm
% Inputs:
%       DRO: the specification of the DRO problem, with fields:
%           stage1: specification of the first-stage problem:
%               c1: the first-stage cost vector
%               Lin: the first-stage inequality constraint matrix (<=)
%               qin: the first-stage inequality constraint rhs
%               Leq: the first-stage equality constraint matrix
%               qeq: the first-stage equality constraint rhs
%               lb: the first-stage variable lower bounds
%               ub: the first-stage variable upper bounds
%           stage2: specification of the second-stage dual problem:
%               V: the matrix linking the first-stage decision and the
%               second-stage objective
%               W: the matrix linking the uncertain parameters and the
%               second-stage objective
%               b: the constant vector in the second-stage objective
%               Lin: the second-stage inequality constraint matrix (<=) 
%               qin: the second-stage inequality constraint rhs
%               Leq: the second-stage equality constraint matrix
%               qeq: the second-stage equality constraint rhs
%               lb: the second-stage variable lower bounds
%               ub: the second-stage variable upper bounds
%       knots: the knots in each dimension (all knots are included)
%       expval: the expected values of the CPWA basis functions in each
%       dimension
%       coup_init: the initial coupling
%           x: matrix of knot indices where each row represents an atom
%           lambda: matrix containing second-stage dual variables, each row
%           represents an atom
%       options: struct with the following fields:
%           tolerance: numerical tolerance (default is 1e-4)
%           reduce_cuts: struct containing options related to reducing cuts
%               thres: constraints (cuts) with slackness above the
%               threshold will be removed (default is inf)
%               freq: frequency of reducing cuts (default is 20)
%               max_iter: no cut will be performed beyond this iteration
%               (default is inf)
%           display: whether to display output (default is true)
%           log_file: path to the log file where the outputs will be
%           written (default is '')
%           LP_params: additional parameters for the LP solver
%           global_params: additional parameters for the global
%           optimization solver
%           rescue: struct containing parameters for saving and loading
%           save files 
%               save_file: path to save the progress (default is '')
%               save_interval: interval to save progress (default is 100)
%               load_file: path to load a save file before running the
%               algorithm (default is '')
% Outputs:
%       LSIP_UB: the computed upper bound for the LSIP problem
%       LSIP_LB: the computed lower bound for the LSIP problem
%       coup_meas: the optimized coupling
%           x: matrix containing atoms of x represented by knot indices
%           lambda: matrix containing atoms of lambda
%           probs: corresponding probabilities
%       stage1_decision: the optimized first-stage decision
%       stage1_multipliers: the optimized multipliers corresponding to the
%       first-stage problem with fields in and eq
%       dual_func: struct representing the dual optimizer
%           y_cell: cell array containing the coefficients of interpolation
%           functions for each measure (the constant intercept is absorbed
%           into the dual function in the first dimension)
%       output: struct containing additional outputs
%           total_time: total time spent in the algorithm
%           iter: number of iterations
%           LP_time: time spent solving LP problems
%           global_time: time spent solving global optimization problems
%       options: return the options in order to retrive some default values

% start the timer
total_timer = tic;

% set the default options
if ~exist('options', 'var') || isempty(options)
    options = struct;
end

if ~isfield(options, 'tolerance') || isempty(options.tolerance)
    options.tolerance = 1e-4;
end

if ~isfield(options, 'reduce_cuts') || isempty(options.reduce_cuts)
    options.reduce_cuts = struct;
end

if ~isfield(options.reduce_cuts, 'thres') ...
        || isempty(options.reduce_cuts.thres)
    options.reduce_cuts.thres = inf;
end

if ~isfield(options.reduce_cuts, 'freq') ...
        || isempty(options.reduce_cuts.freq)
    options.reduce_cuts.freq = 20;
end

if ~isfield(options.reduce_cuts, 'max_iter') ...
        || isempty(options.reduce_cuts.max_iter)
    options.reduce_cuts.max_iter = inf;
end

if ~isfield(options, 'display') || isempty(options.display)
    options.display = true;
end

if ~isfield(options, 'log_file') || isempty(options.log_file)
    options.log_file = '';
end

if ~isfield(options, 'rescue') || isempty(options.rescue)
    options.rescue = struct;
end

if ~isfield(options.rescue, 'save_file') ...
        || isempty(options.rescue.save_file)
    options.rescue.save_file = '';
end

if ~isfield(options.rescue, 'save_interval') ...
        || isempty(options.rescue.save_interval)
    options.rescue.save_interval = 100;
end

if ~isfield(options.rescue, 'load_file') ...
        || isempty(options.rescue.load_file)
    options.rescue.load_file = '';
end

% dimensionality
K1 = length(DRO.stage1.c1);
K2ast = size(DRO.stage2.V, 1);
N = length(knots);

assert(length(expval) == N, 'expected values invalid');


stage1_in_no = 0;
stage1_eq_no = 0;

% check the DRO specifications
if isfield(DRO.stage1, 'Lin') && isfield(DRO.stage1, 'qin')
    assert(size(DRO.stage1.Lin, 2) == K1);
    assert(size(DRO.stage1.Lin, 1) == length(DRO.stage1.qin));
    stage1_in_no = length(DRO.stage1.qin);
end

if isfield(DRO.stage1, 'Leq') && isfield(DRO.stage1, 'qeq')
    assert(size(DRO.stage1.Leq, 2) == K1);
    assert(size(DRO.stage1.Leq, 1) == length(DRO.stage1.qeq));
    stage1_eq_no = length(DRO.stage1.qeq);
end

if isfield(DRO.stage1, 'lb')
    assert(length(DRO.stage1.lb) == K1);
    stage1_in_no = stage1_in_no + sum(~isinf(DRO.stage1.lb));
end

if isfield(DRO.stage1, 'ub')
    assert(length(DRO.stage1.ub) == K1);
    stage1_in_no = stage1_in_no + sum(~isinf(DRO.stage1.ub));
end

assert(size(DRO.stage2.V, 1) == K2ast);
assert(size(DRO.stage2.V, 2) == K1);
assert(size(DRO.stage2.W, 1) == K2ast);
assert(size(DRO.stage2.W, 2) == N);
assert(length(DRO.stage2.b) == K2ast);

if isfield(DRO.stage2, 'Lin') && isfield(DRO.stage2, 'qin')
    assert(size(DRO.stage2.Lin, 2) == K2ast);
    assert(size(DRO.stage2.Lin, 1) == length(DRO.stage2.qin));
end

if isfield(DRO.stage2, 'Leq') && isfield(DRO.stage2, 'qeq')
    assert(size(DRO.stage2.Leq, 2) == K2ast);
    assert(size(DRO.stage2.Leq, 1) == length(DRO.stage2.qeq));
end

assert(length(DRO.stage2.lb) == K2ast);
assert(length(DRO.stage2.ub) == K2ast);

% number of knots in each dimension
knot_no_list = zeros(N, 1);

% the list used to remove the first know from dimensions 2 to N
keep_list_cell = cell(N, 1);

% initialize the stage1 multipliers
stage1_multipliers = struct;
stage1_multipliers.in = zeros(stage1_in_no, 1);
stage1_multipliers.eq = zeros(stage1_eq_no, 1);

% the objective vectors where the first basis function is removed for
% identification
v_cell = cell(N, 1);

% initialize the dual optimizer
dual_func = struct;
dual_func.y_cell = cell(N, 1);

for i = 1:N
    knot_no_list(i) = length(knots{i});
    keep_list_cell{i} = true(knot_no_list(i), 1);
    
    if i == 1
        % the constant is absorbed into the first dimension
        v_cell{i} = expval{i};
        dual_func.y_cell{i} = zeros(knot_no_list(i), 1);
    else
        % the first basis function is removed for identification
        keep_list_cell{i}(1) = false;
        v_cell{i} = expval{i}(2:end);
        dual_func.y_cell{i} = zeros(knot_no_list(i) - 1, 1);
    end
    
    assert(size(expval{i}, 1) == knot_no_list(i), ...
        'knots and expected values mismatch');
end

% some data structures to simplify computation later
knot_vec = vertcat(knots{:});
keep_list = vertcat(keep_list_cell{:});

% the offset of knots in each dimension
knot_offset = [0; cumsum(knot_no_list)];
knot_offset = knot_offset(1:end - 1);

% the offset of points in each dimension
var_offset = [0; cumsum(knot_no_list(1:end - 1) - 1) + 1] + K1;

% total number of decision variables in the LP problem including y0
LP_var_len = K1 + sum(knot_no_list - 1) + 1;

% assemble the Gurobi model for the linear programming problem
% decision variable: [a; y] where a is the first-stage decision and y
% represents the coefficients of the interpolation basis functions
LP_model = struct;
LP_model.modelsense = 'min';
LP_model.obj = [DRO.stage1.c1; vertcat(v_cell{:})];
LP_model.lb = -inf(LP_var_len, 1);
LP_model.ub = inf(LP_var_len, 1);

% generate the constraints for the first-stage decision
stage1_A = sparse(0, K1);
stage1_rhs = zeros(0, 1);

if isfield(DRO.stage1, 'Lin') && isfield(DRO.stage1, 'qin')
    % append the inequality constraints
    stage1_A = [stage1_A; sparse(DRO.stage1.Lin)];
    stage1_rhs = [stage1_rhs; DRO.stage1.qin];
end

if isfield(DRO.stage1, 'lb')
    % turn lower bounds into inequality constraints
    lb_constrained = ~isinf(DRO.stage1.lb);
    stage1_A = [stage1_A; sparse(1:sum(lb_constrained), ...
        find(lb_constrained), -1, sum(lb_constrained), K1)];
    stage1_rhs = [stage1_rhs; -DRO.stage1.lb(lb_constrained)];
end

if isfield(DRO.stage1, 'ub')
    % turn upper bounds into inequality constraints
    ub_constrained = ~isinf(DRO.stage1.ub);
    stage1_A = [stage1_A; sparse(1:sum(ub_constrained), ...
        find(ub_constrained), 1, sum(ub_constrained), K1)];
    stage1_rhs = [stage1_rhs; DRO.stage1.ub(ub_constrained)];
end

if isfield(DRO.stage1, 'Leq') && isfield(DRO.stage1, 'qeq')
    % append the equality constraints (equality constraints are placed
    % after inequality constraints)
    stage1_A = [stage1_A; sparse(DRO.stage1.Leq)];
    stage1_rhs = [stage1_rhs; DRO.stage1.qeq];
end

stage1_sense = [repmat('<', stage1_in_no, 1); 
    repmat('=', stage1_eq_no, 1)];

% generate the initial constraints coming from the initial coupling
[A_init, rhs_init] = DRO_feascons(coup_init.x, coup_init.lambda, ...
    knot_vec, knot_offset, keep_list, DRO);
sense_init = repmat('>', length(rhs_init), 1);
LP_model.A = [stage1_A, sparse(length(stage1_rhs), LP_var_len - K1); 
    A_init];
LP_model.rhs = [stage1_rhs; rhs_init];
LP_model.sense = [stage1_sense; sense_init];


% parameters of the LP solver
LP_params = struct;

% disable output by default
LP_params.OutputFlag = 0;

% since the gurobi MATLAB interface does not support callback, set a
% 30-minute time limit on the LP solver; if the solver does not converge
% when the time limit is hit, it assumes that numerical issues have
% occurred and restarts the solution process with NumericFocus = 3 and
% LPWarmStart = 2 and TimeLimit = inf
LP_params.TimeLimit = 1800;

% set the additional parameters for the LP solver
if isfield(options, 'LP_params') && ~isempty(options.LP_params)
    LP_params_fields = fieldnames(options.LP_params);
    LP_params_values = struct2cell(options.LP_params);
    
    for fid = 1:length(LP_params_fields)
        LP_params.(LP_params_fields{fid}) = LP_params_values{fid};
    end
end

% parameters of the global (MILP) solver
global_params = struct;
global_params.OutputFlag = 0;
global_params.IntFeasTol = 1e-6;
global_params.FeasibilityTol = 1e-8;
global_params.OptimalityTol = 1e-8;
global_params.PoolSolutions = 100;
global_params.PoolGap = 0.8;
global_params.NodefileStart = 2;
global_params.BestBdStop = 1e-6;
global_params.BestObjStop = -inf;
global_params.MIPGap = 1e-4;
global_params.MIPGapAbs = 1e-10;

% set the additional parameters for the global (MILP) solver
if isfield(options, 'global_params') ...
        && ~isempty(options.global_params)
    global_params_fields = fieldnames(options.global_params);
    global_params_values = struct2cell(options.global_params);
    
    for fid = 1:length(global_params_fields)
        global_params.(global_params_fields{fid}) ...
            = global_params_values{fid};
    end
end

% compute the truncation bounds for the conjugate functions
LP_conjbd_model = DRO_stage2_gurobi(DRO);

conjbd = zeros(N, 2);

for i = 1:N
    for j = 1:2
        LP_conjbd_model.obj = (j * 2 - 3) * full(DRO.stage2.W(:, i));

        LP_conjbd_output = gurobi(LP_conjbd_model, LP_params);
        conjbd(i, j) = (j * 2 - 3) * LP_conjbd_output.objval;
    end
end

% initialize the lower and upper bounds
LSIP_UB = inf;
LSIP_LB = -inf;
iter = 0;

% initialize the set of constraints
x_atoms_agg = coup_init.x;
lambda_atoms_agg = coup_init.lambda;

% the basis used for warm start
cbasis = [];
vbasis = [];

% some statistics
LP_time = 0;
global_time = 0;
total_time = 0;

% open the log file
if ~isempty(options.log_file)
    log_file = fopen(options.log_file, 'a');

    if log_file < 0
        error('cannot open log file');
    end
end

if ~isempty(options.rescue.load_file)
    % load all states from the load file
    load(options.rescue.load_file);
end

% loop until the gap between lower and upper bounds is below the tolerance
while true
    % solve the LP problem
    LP_timer = tic;
    LP_trial_num = 0;
    while true
        LP_params_runtime = LP_params;

        if LP_trial_num == 1
            % if the LP solver has failed once (reaching the time limit
            % without converging), retry with high numeric focus and
            % presolve (with warm-start)
            LP_params_runtime.TimeLimit = inf;
            LP_params_runtime.NumericFocus = 3;
            LP_params_runtime.LPWarmStart = 2;
        end

        LP_result = gurobi(LP_model, LP_params_runtime);

        if strcmp(LP_result.status, 'OPTIMAL')
            break;
        else
            LP_trial_num = LP_trial_num + 1;
        end

        if LP_trial_num >= 2
            % if the LP solver fails after the second trial, report error
            error('unexpected error while solving LP');
        end
    end
    LP_time = LP_time + toc(LP_timer);
    
    if isfield(LP_result, 'vbasis') && ~isempty(LP_result.vbasis) ...
            && isfield(LP_result, 'cbasis') && ~isempty(LP_result.cbasis)
        vbasis = LP_result.vbasis;
        cbasis = LP_result.cbasis;
    end
    
    % store the LP dual optimizer
    LP_dual = LP_result.pi;
    
    % update the lower bound
    LSIP_LB = LP_result.objval;

    % store the first-stage decision and multipliers
    stage1_decision = LP_result.x(1:K1);
    stage1_multipliers.in = LP_dual(1:stage1_in_no);
    stage1_multipliers.eq = LP_dual(stage1_in_no + (1:stage1_eq_no));
    
    % compute the value of the CPWA functions (including the first knot in
    % each dimension after the first)
    val_cell = cell(N, 1);
    for i = 1:N
        if i == 1
            val_cell{i} = LP_result.x(var_offset(i) + (1:knot_no_list(i)));
        else
            val_cell{i} = [0; LP_result.x(var_offset(i) ...
                + (1:knot_no_list(i) - 1));];
        end
    end
    
    % solve the global optimization problem

    % prepare the inputs for the global optimization problem
    env_knots = cell(N, 1);
    env_vals = cell(N, 1);
    env_vals_old = cell(N, 1);
    concave_knots = cell(N, 1);
    concave_val = cell(N, 1);

    % the objective value of a feasible dual (LSIP) solution
    LSIP_UB_cand = stage1_decision' * DRO.stage1.c1;
    
    for i = 1:N
        % take the convex envelope function
        [env_knots{i}, env_vals{i}, env_vals_old{i}] ...
            = CPWA_1d_convexenv(knots{i}, val_cell{i});

        % store the basis representation of the convex envelope function
        if i == 1
            dual_func.y_cell{i} = env_vals_old{i};
        else
            dual_func.y_cell{i} = env_vals_old{i}(2:end);
        end
        LSIP_UB_cand = LSIP_UB_cand + dual_func.y_cell{i}' * v_cell{i};
        
        % compute the convex conjugate functions
        [concave_knots{i}, concave_val{i}] = CPWA_1d_conjugate( ...
            env_knots{i}, env_vals{i}, conjbd(i, :)');
        concave_val{i} = -concave_val{i};
    end
    
    global_model = CPWA_sep_concavemin_MILP_gurobi(0, ...
        -DRO.stage2.V * stage1_decision - DRO.stage2.b, ...
        DRO.stage2.W', concave_knots, concave_val, DRO.stage2);
    global_timer = tic;
    global_result = gurobi(global_model, global_params);
    global_time = global_time + toc(global_timer);

    % absorb the added constant into the dual function in the first
    % dimension
    dual_func.y_cell{1} = dual_func.y_cell{1} - global_result.objbound;

    LSIP_UB_cand = LSIP_UB_cand - global_result.objbound;
    LSIP_UB = min(LSIP_UB, LSIP_UB_cand);

    % display output
    if options.display
        fprintf('iter = %6d, LB = %10.6f, UB = %10.6f\n', iter, ...
            LSIP_LB, LSIP_UB);
    end

    % write log
    if ~isempty(options.log_file)
        fprintf(log_file, 'iter = %6d, LB = %10.6f, UB = %10.6f\n', ...
            iter, LSIP_LB, LSIP_UB);
    end

    % check the termination criterion here
    if LSIP_UB - LSIP_LB <= options.tolerance
        break;
    end
    
    % reduce cuts
    if ~isinf(options.reduce_cuts.thres) ...
            && iter > 0 && iter <= options.reduce_cuts.max_iter ...
            && mod(iter, options.reduce_cuts.freq) == 0
        % the list of constraints to be kept
        constr_keep_list = LP_result.slack >= -options.reduce_cuts.thres;

        % always keep the constraints the initial constraints including the
        % constraints on the first-stage decision
        constr_keep_list(1:(length(rhs_init) + length(stage1_rhs))) = true;
        cut_keep_list = constr_keep_list(length(stage1_rhs) + 1:end);
        
        % update all variables
        x_atoms_agg = x_atoms_agg(cut_keep_list, :);
        lambda_atoms_agg = lambda_atoms_agg(cut_keep_list, :);
        LP_model.A = LP_model.A(constr_keep_list, :);
        LP_model.rhs = LP_model.rhs(constr_keep_list);
        LP_model.sense = LP_model.sense(constr_keep_list);
        
        if ~isempty(cbasis)
            cbasis = cbasis(constr_keep_list);
        end
    end
    
    % get a set of approximate optimizers of the global optimization
    % problem
    if isfield(global_result, 'pool')
        lambda_cut = horzcat(global_result.pool.xn);
        lambda_cut = lambda_cut(1:K2ast, [global_result.pool.objval] < 0)';
    else
        lambda_cut = global_result.x(1:K2ast)';
    end

    % clip the lambda atoms to the bounds to resolve small numerical issues
    lambda_cut = min(max(lambda_cut, DRO.stage2.lb'), DRO.stage2.ub');

    trans_cut = (DRO.stage2.W' * lambda_cut')';
        
    % clip the transformed vector at the lower and upper bounds to
    % resolve small numerical issues
    trans_cut = min(max(trans_cut, conjbd(:, 1)'), conjbd(:, 2)');

    x_cut = zeros(size(trans_cut, 1), N);

    % retrieve the x knot indices
    for i = 1:N
        [~, x_cut(:, i)] = min(env_vals_old{i}' - trans_cut(:, i) ...
            .* knots{i}', [], 2);
    end

    [A_new, rhs_new] = DRO_feascons(x_cut, lambda_cut, ...
        knot_vec, knot_offset, keep_list, DRO);
    sense_new = repmat('>', length(rhs_new), 1);

    % generate new constraints
    LP_model.A = [LP_model.A; A_new];
    LP_model.rhs = [LP_model.rhs; rhs_new];
    LP_model.sense = [LP_model.sense; sense_new];
    x_atoms_agg = [x_atoms_agg; x_cut]; %#ok<AGROW> 
    lambda_atoms_agg = [lambda_atoms_agg; lambda_cut]; %#ok<AGROW> 
    
    if ~isempty(vbasis) && ~isempty(cbasis)
        LP_model.vbasis = vbasis;
        LP_model.cbasis = [cbasis; zeros(length(rhs_new), 1)];
    else
        LP_model = rmfield(LP_model, 'vbasis');
        LP_model = rmfield(LP_model, 'cbasis');
    end
    
    iter = iter + 1;

    % overwrite manual_save to true to save states while debugging
    manual_save = false;

    % save states
    if ~isempty(options.rescue.save_file) ...
            && (mod(iter, options.rescue.save_interval) == 0 ...
            || manual_save)
        save(options.rescue.save_file, 'A_init', 'A_new', ...
            'cbasis', 'concave_knots', 'concave_val', 'conjbd', ...
            'coup_init', 'DRO', 'dual_func', 'env_knots', 'env_vals', ...
            'env_vals_old', 'expval', 'global_model', 'global_result', ...
            'global_time', 'iter', 'K1', 'K2ast', 'keep_list', ...
            'keep_list_cell', 'knot_no_list', 'knot_offset', ...
            'knot_vec', 'knots', 'lambda_atoms_agg', 'lambda_cut', ...
            'LP_conjbd_model', 'LP_conjbd_output', ...
            'LP_dual', 'LP_model', 'LP_result', 'LP_time', ...
            'LP_var_len', 'LSIP_LB', 'LSIP_UB', 'LSIP_UB_cand', ...
            'N', 'rhs_init', 'rhs_new', 'sense_init', 'sense_new', ...
            'stage1_A', 'stage1_decision', 'stage1_eq_no', ...
            'stage1_in_no', 'stage1_multipliers', 'stage1_rhs', ...
            'stage1_sense', 'total_time', 'trans_cut', ...
            'v_cell', 'val_cell', 'var_offset', ...
            'vbasis', 'x_atoms_agg', 'x_cut', ...
            '-v7.3');
    end
end

% the approxiamte optimal coupling
LP_dual_semiinf = LP_dual(stage1_in_no + stage1_eq_no + 1:end);
% remove atoms with very small probabilities
probs_nonzero = LP_dual_semiinf > eps;
coup_meas = struct;
coup_meas.x = x_atoms_agg(probs_nonzero, :);
coup_meas.lambda = lambda_atoms_agg(probs_nonzero, :);
coup_meas.probs = LP_dual_semiinf(probs_nonzero);
% normalize the probabilities in case there is some numerical inaccuracies
coup_meas.probs = coup_meas.probs / sum(coup_meas.probs);


% stop the timer
total_time = total_time + toc(total_timer);

% prepare additional outputs
output = struct;

output.total_time = total_time;
output.iter = iter;
output.LP_time = LP_time;
output.global_time = global_time;

% close the log file
if ~isempty(options.log_file)
    fclose(log_file);   
end

end