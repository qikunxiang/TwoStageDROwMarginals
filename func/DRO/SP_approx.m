function [objval, stage1_decision, output] = SP_approx(SP, samp, ...
    weight, options)
% Approximately solve a stochastic programming problem by the sample
% average approximation
% Inputs: 
%       SP: the specification of the SP problem, with fields:
%           stage1: specification of the first-stage problem:
%               c1: the first-stage cost vector
%               Lin: the first-stage inequality constraint matrix (<=)
%               qin: the first-stage inequality constraint rhs
%               Leq: the first-stage equality constraint matrix
%               qeq: the first-stage equality constraint rhs
%               lb: the first-stage variable lower bounds
%               ub: the first-stage variable upper bounds
%           stage2: specification of the second-stage dual problem:
%               c2: the second-stage cost vector
%               Ain: the second-stage inequality contraint matrix (<=)
%               Vin: the matrix linking the first-stage decision and the
%               rhs of the second-stage inequality constraints
%               Win: the matrix linking the uncertain parameters and the
%               rhs of the second-stage inequality constraints
%               bin: the constant vector in the rhs of the second-stage 
%               inequality constraints
%               Aeq: the second-stage equality contraint matrix
%               Veq: the matrix linking the first-stage decision and the
%               rhs of the second-stage equality constraints
%               Weq: the matrix linking the uncertain parameters and the
%               rhs of the second-stage equality constraints
%               beq: the constant vector in the rhs of the second-stage 
%               equality constraints
%               lb: the second-stage variable lower bounds
%               ub: the second-stage variable upper bounds
%       samp: matrix containing samples of the uncertain parameters in rows
%       weight: probability weights associated with the samples
%       options: struct with the following fields:
%           log_file: path to the log file where the outputs will be
%           written (default is '')
%           LP_params: additional parameters for the LP solver
% Outputs:
%       objval: the minimized objective
%       stage1_decision: the optimized first-stage decision
%       output: struct containing additional outputs
%           total_time: total time spent in the algorithm
%           LP_time: time spent solving LP problems

% start the timer
total_timer = tic;

% set the default options
if ~exist('options', 'var') || isempty(options)
    options = struct;
end

if ~isfield(options, 'log_file') || isempty(options.log_file)
    options.log_file = '';
end

% dimensionality
K1 = length(SP.stage1.c1);
K2 = length(SP.stage2.c2);
N = size(samp, 2);
samp_no = size(samp, 1);

assert(length(weight) == samp_no, 'sample weights invalid');

% check the DRO specifications
if isfield(SP.stage1, 'Lin') && isfield(SP.stage1, 'qin')
    assert(size(SP.stage1.Lin, 2) == K1);
    assert(size(SP.stage1.Lin, 1) == length(SP.stage1.qin));
    stage1_in_no = length(SP.stage1.qin);
else
    stage1_in_no = 0;
end

if isfield(SP.stage1, 'Leq') && isfield(SP.stage1, 'qeq')
    assert(size(SP.stage1.Leq, 2) == K1);
    assert(size(SP.stage1.Leq, 1) == length(SP.stage1.qeq));
    stage1_eq_no = length(SP.stage1.qeq);
else
    stage1_eq_no = 0;
end

if isfield(SP.stage1, 'lb')
    assert(length(SP.stage1.lb) == K1);
end

if isfield(SP.stage1, 'ub')
    assert(length(SP.stage1.ub) == K1);
end

if isfield(SP.stage2, 'Ain') && isfield(SP.stage2, 'Vin') ...
        && isfield(SP.stage2, 'Win') && isfield(SP.stage2, 'bin')
    stage2_in_no = size(SP.stage2.Ain, 1);
    assert(size(SP.stage2.Ain, 2) == K2);
    assert(size(SP.stage2.Vin, 1) == stage2_in_no);
    assert(size(SP.stage2.Vin, 2) == K1);
    assert(size(SP.stage2.Win, 1) == stage2_in_no);
    assert(size(SP.stage2.Win, 2) == N);
    assert(length(SP.stage2.bin) == stage2_in_no);
else
    stage2_in_no = 0;
end

if isfield(SP.stage2, 'Aeq') && isfield(SP.stage2, 'Veq') ...
        && isfield(SP.stage2, 'Weq') && isfield(SP.stage2, 'beq')
    stage2_eq_no = size(SP.stage2.Aeq, 1);
    assert(size(SP.stage2.Aeq, 2) == K2);
    assert(size(SP.stage2.Veq, 1) == stage2_eq_no);
    assert(size(SP.stage2.Veq, 2) == K1);
    assert(size(SP.stage2.Weq, 1) == stage2_eq_no);
    assert(size(SP.stage2.Weq, 2) == N);
    assert(length(SP.stage2.beq) == stage2_eq_no);
else
    stage2_eq_no = 0;
end

if isfield(SP.stage2, 'lb')
    assert(length(SP.stage2.lb) == K2);
end

if isfield(SP.stage2, 'ub')
    assert(length(SP.stage2.ub) == K2);
end

% build the LP model for stochastic programming with weighted samples
LP_model = struct;
LP_model.modelsense = 'min';
LP_model.obj = [SP.stage1.c1; reshape(SP.stage2.c2 .* weight', [], 1)];

% first-stage inequality constraints
if stage1_in_no > 0
    A_stage1_in = [sparse(SP.stage1.Lin), sparse(stage1_in_no, ...
        K2 * samp_no)];
    rhs_stage1_in = SP.stage1.qin;
    sense_stage1_in = repmat('<', stage1_in_no, 1);
else
    A_stage1_in = sparse(0, K1 + K2 * samp_no);
    rhs_stage1_in = zeros(0, 1);
    sense_stage1_in = [];
end

% first-stage equality constraints
if stage1_eq_no > 0
    A_stage1_eq = [sparse(SP.stage1.Leq), sparse(stage1_eq_no, ...
        K2 * samp_no)];
    rhs_stage1_eq = SP.stage1.qeq;
    sense_stage1_eq = repmat('=', stage1_eq_no, 1);
else
    A_stage1_eq = sparse(0, K1 + K2 * samp_no);
    rhs_stage1_eq = zeros(0, 1);
    sense_stage1_eq = [];
end

% second-stage inequality constraints
if stage2_in_no > 0
    A_stage2_in = [repmat(-sparse(SP.stage2.Vin), samp_no, 1), ...
        kron(speye(samp_no), sparse(SP.stage2.Ain))];
    rhs_stage2_in = full(reshape(SP.stage2.Win * samp' ...
        + SP.stage2.bin, [], 1));
    sense_stage2_in = repmat('<', samp_no * stage2_in_no, 1);
else
    A_stage2_in = sparse(0, K1 + K2 * samp_no);
    rhs_stage2_in = zeros(0, 1);
    sense_stage2_in = [];
end

% second-stage equality constraints
if stage2_eq_no > 0
    A_stage2_eq = [repmat(-sparse(SP.stage2.Veq), samp_no, 1), ...
        kron(speye(samp_no), sparse(SP.stage2.Aeq))];
    rhs_stage2_eq = full(reshape(SP.stage2.Weq * samp' ...
        + SP.stage2.beq, [], 1));
    sense_stage2_eq = repmat('=', samp_no * stage2_eq_no, 1);
else
    A_stage2_eq = sparse(0, K1 + K2 * samp_no);
    rhs_stage2_eq = zeros(0, 1);
    sense_stage2_eq = [];
end

LP_model.A = [A_stage1_in; A_stage1_eq; A_stage2_in; A_stage2_eq];
LP_model.rhs = [rhs_stage1_in; rhs_stage1_eq; 
    rhs_stage2_in; rhs_stage2_eq];
LP_model.sense = [sense_stage1_in; sense_stage1_eq; 
    sense_stage2_in; sense_stage2_eq];

% bounds
if isfield(SP.stage1, 'lb')
    lb_stage1 = SP.stage1.lb;
else
    lb_stage1 = -inf(K1, 1);
end

if isfield(SP.stage1, 'ub')
    ub_stage1 = SP.stage1.ub;
else
    ub_stage1 = inf(K1, 1);
end

if isfield(SP.stage2, 'lb')
    lb_stage2 = SP.stage2.lb;
else
    lb_stage2 = -inf(K2, 1);
end

if isfield(SP.stage2, 'ub')
    ub_stage2 = SP.stage2.ub;
else
    ub_stage2 = inf(K2, 1);
end

LP_model.lb = [lb_stage1; repmat(lb_stage2, samp_no, 1)];
LP_model.ub = [ub_stage1; repmat(ub_stage2, samp_no, 1)];

% parameters of the LP solver
LP_params = struct;

% disable output by default
LP_params.OutputFlag = 0;

% set the additional parameters for the LP solver
if isfield(options, 'LP_params') && ~isempty(options.LP_params)
    LP_params_fields = fieldnames(options.LP_params);
    LP_params_values = struct2cell(options.LP_params);
    
    for fid = 1:length(LP_params_fields)
        LP_params.(LP_params_fields{fid}) = LP_params_values{fid};
    end
end

% open the log file
if ~isempty(options.log_file)
    log_file = fopen(options.log_file, 'a');
    
    if log_file < 0
        error('cannot open log file');
    end
end

% solve the LP
LP_timer = tic;
LP_output = gurobi(LP_model, LP_params);
LP_time = toc(LP_timer);

if ~strcmp(LP_output.status, 'OPTIMAL')
    error('LP solver did not find an optimal solution');
end

objval = LP_output.objval;
stage1_decision = LP_output.x(1:K1);

% stop the timer
total_time = toc(total_timer);

% prepare additional outputs
output = struct;

output.total_time = total_time;
output.LP_time = LP_time;

% close the log file
if ~isempty(options.log_file)
    fclose(log_file);   
end

end

