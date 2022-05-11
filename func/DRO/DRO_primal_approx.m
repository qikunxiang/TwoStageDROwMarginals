function val_rep = DRO_primal_approx(DRO, stage1_multipliers, ...
    coup_meas, marginv, samp_no, rep_no, batch_size, rand_stream, ...
    samp_mode, run_in_parallel)
% Compute a lower bound on the two-stage DRO problem with marginal 
% constraints via Monte Carlo sampling from a coupling when the first-stage 
% decision is not fixed
% Inputs: 
%       DRO: the specification of the DRO problem
%       stage1_multipliers: the optimized multipliers corresponding to the
%       first-stage problem with fields in and eq
%       coup_meas: the given coupling
%       marginv: the inverse cdf of the marginal distributions represented
%       using a single function handle (i, u) -> marginv(i, u)
%       samp_no: the number of samples used in the Monte Carlo integration
%       rep_no: the number of repetitions to approximate the distribution 
%       of the Monte Carlo integration
%       batch_size: the maximum number of samples generated in each batch,
%       used to prevent memory surges (default is 1e6)
%       rand_stream: RandStream object to guarantee reproducibility
%       samp_mode: the copula used for generating uniform random variables,
%       either 'comonotone' or 'independent' (default is 'comonotone')
%       run_in_parallel: a boolean parameter indicating whether to compute 
%       the repetitions in parallel (default is false)
% Outputs: 
%       val_rep: a list of lower bounds from repeated trials of
%       Monte Carlo integration

if ~exist('batch_size', 'var') || isempty(batch_size)
    batch_size = 1e6;
end

if ~exist('samp_mode', 'var') || isempty(samp_mode)
    samp_mode = 'comonotone';
end

if ~exist('run_in_parallel', 'var') || isempty(run_in_parallel)
    run_in_parallel = false;
end

K1 = length(DRO.stage1.c1);
K2ast = size(DRO.stage2.V, 1);
N = size(DRO.stage2.W, 2);

% check the DRO specifications
stage1_in_no = 0;
stage1_in_rhs = zeros(0, 1);
stage1_eq_no = 0;
stage1_eq_rhs = zeros(0, 1);

if isfield(DRO.stage1, 'Lin') && isfield(DRO.stage1, 'qin')
    assert(size(DRO.stage1.Lin, 2) == K1);
    assert(size(DRO.stage1.Lin, 1) == length(DRO.stage1.qin));
    stage1_in_no = length(DRO.stage1.qin);
    stage1_in_rhs = [stage1_in_rhs; DRO.stage1.qin];
end

if isfield(DRO.stage1, 'Leq') && isfield(DRO.stage1, 'qeq')
    assert(size(DRO.stage1.Leq, 2) == K1);
    assert(size(DRO.stage1.Leq, 1) == length(DRO.stage1.qeq));
    stage1_eq_no = length(DRO.stage1.qeq);
    stage1_eq_rhs = [stage1_eq_rhs; DRO.stage1.qeq];
end

if isfield(DRO.stage1, 'lb')
    assert(length(DRO.stage1.lb) == K1);
    stage1_lb_active = ~isinf(DRO.stage1.lb);
    stage1_in_no = stage1_in_no + sum(stage1_lb_active);
    stage1_in_rhs = [stage1_in_rhs; -DRO.stage1.lb(stage1_lb_active)];
end

if isfield(DRO.stage1, 'ub')
    assert(length(DRO.stage1.ub) == K1);
    stage1_ub_active = ~isinf(DRO.stage1.ub);
    stage1_in_no = stage1_in_no + sum(stage1_ub_active);
    stage1_in_rhs = [stage1_in_rhs; DRO.stage1.ub(stage1_ub_active)];
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

% check the inputs
atom_no = length(coup_meas.probs);
assert(size(coup_meas.x, 2) == N && size(coup_meas.x, 1) == atom_no ...
    && size(coup_meas.lambda, 2) == K2ast && size(coup_meas.lambda, 1) ...
    == atom_no, 'input measure mis-specified'); 

% start the computation
val_rep = zeros(rep_no, 1);

% compute the term in the integral that does not require Monte Carlo
const_part = 0;

if stage1_in_no > 0
    const_part = const_part + stage1_multipliers.in' * stage1_in_rhs;
end

if stage1_eq_no > 0
    const_part = const_part + stage1_multipliers.eq' * stage1_eq_rhs;
end

if run_in_parallel

    rs = parallel.pool.Constant(rand_stream);

    parfor_progress(rep_no);
    parfor rep_id = 1:rep_no
        substream = rs.Value;
        substream.Substream = rep_id;
        
        val_sum = 0;
        samp_no_remain = samp_no;
        
        % cap the number of samples generated at once at the batch size
        while samp_no_remain > 0
            samp_no_batch = min(samp_no_remain, batch_size);
            
            % generate samples for computing the global lower bound
            [copula_samples, atom_list] = discretecopula_samp( ...
                coup_meas.x, coup_meas.probs, samp_no_batch, ...
                substream, samp_mode); %#ok<PFBNS>

            x_samples = zeros(samp_no_batch, N);

            for i = 1:N
                x_samples(:, i) = marginv(i, copula_samples(:, i)); ...
                    %#ok<PFBNS>
            end
            
            lambda_samples = coup_meas.lambda(atom_list, :);
            
            val_samples = sum((DRO.stage2.W * x_samples' ...
                + DRO.stage2.b) .* lambda_samples', 1); %#ok<PFBNS>
            
            % accumulate the samples
            val_sum = val_sum + sum(val_samples);
            
            % subtract the number of samples generated
            samp_no_remain = samp_no_remain - samp_no_batch;
        end
        
        val_rep(rep_id) = val_sum / samp_no + const_part;
        parfor_progress;
    end
    parfor_progress(0);
else
    parfor_progress(rep_no);
    for rep_id = 1:rep_no
        rand_stream.Substream = rep_id;
        
        val_sum = 0;
        samp_no_remain = samp_no;
        
        % cap the number of samples generated at once at the batch size
        while samp_no_remain > 0
            samp_no_batch = min(samp_no_remain, batch_size);
            
            % generate samples for computing the global lower bound
            [copula_samples, atom_list] = discretecopula_samp( ...
                coup_meas.x, coup_meas.probs, samp_no_batch, ...
                rand_stream, samp_mode);

            x_samples = zeros(samp_no_batch, N);

            for i = 1:N
                x_samples(:, i) = marginv(i, copula_samples(:, i));
            end
            
            lambda_samples = coup_meas.lambda(atom_list, :);
            
            val_samples = sum((DRO.stage2.W * x_samples' ...
                + DRO.stage2.b) .* lambda_samples', 1);
            
            % accumulate the samples
            val_sum = val_sum + sum(val_samples);
            
            % subtract the number of samples generated
            samp_no_remain = samp_no_remain - samp_no_batch;
        end
        
        val_rep(rep_id) = val_sum / samp_no + const_part;
        parfor_progress;
    end
    parfor_progress(0);
end

end

