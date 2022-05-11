rng(1000, 'combRecursive');

N = 20;
T = 20;

K1 = N;
K2ast = N;

% problem specification

% the mu and sigma^2 parameters of the mixture of truncated normal
% marginals 
mixtrnorm_trunc = [0; 2];
mixtrnorm_mu = [0.1; 0.5; 1];
mixtrnorm_sig2 = [0.25; 0.04; 0.01];
mixtrnorm_w = [0.7; 0.2; 0.1];

DRO = DRO_task_scheduling(N, T);

% approximation scheme
% the tolerance value used when solving LSIP
LSIP_tolerance = 1e-3;

% approximation scheme
knot_no_list = round(exp(linspace(log(5), log(100), 8)))';
step_no = length(knot_no_list);

knots_cell = cell(step_no, 1);
expval_cell = cell(step_no, 1);
bd_cell = cell(step_no, 1);
err_bound_list = zeros(step_no, 1);
coef_num_cell = cell(step_no, 1);

% first generate the maximum number of knots
[~, knots_hist, bd_hist] = mixtruncnorm_momentset_greedy( ...
    mixtrnorm_mu, mixtrnorm_sig2, mixtrnorm_w, mixtrnorm_trunc, ...
    knot_no_list(end));

for step_id = 1:step_no
    knots = cell(N, 1);
    expval = cell(N, 1);
    bd_list = zeros(N, 1);
    coef_num_list = zeros(N, 1);
    
    for i = 1:N
        % retrive a subset of knots
        knots{i} = sort(knots_hist(1:knot_no_list(step_id)), ...
            'ascend');

        % the corresponding upper bound on the Wasserstein-1 radius
        bd_list(i) = bd_hist(knot_no_list(step_id));
        
        % compute the corresponding integrals
        expval{i} = mixtruncnorm_momentset(mixtrnorm_mu, ...
            mixtrnorm_sig2, mixtrnorm_w, knots{i});

        % the number of coefficients (with the constant intercept)
        coef_num_list(i) = length(knots{i});
    end
    
    knots_cell{step_id} = knots;
    expval_cell{step_id} = expval;
    err_bound_list(step_id) = LSIP_tolerance + sum(bd_list) * N;
    coef_num_cell{step_id} = coef_num_list;
end

marg1 = cell(N, 2);

for i = 1:N
    % generate a discrete measure in the moment set
    [marg1{i, 1}, marg1{i, 2}] = momentset1d_measinit_bounded( ...
        (1:length(knots_cell{1}{i}))', expval_cell{1}{i});
end

% compute the comonotone coupling as the initial measure
[joint_atoms, joint_probs] = comonotone_coupling(marg1);

coup_init1 = struct;
coup_init1.probs = joint_probs;
coup_init1.x = joint_atoms;

save('exp/DRO_scheduling/exp_inputs.mat', ...
    'N', 'K1', 'K2ast', 'DRO', ...
    'mixtrnorm_mu', 'mixtrnorm_sig2', 'mixtrnorm_w', ...
    'mixtrnorm_trunc', 'LSIP_tolerance', 'err_bound_list', ...
    'knot_no_list', 'step_no', 'knots_cell', 'expval_cell', ...
    'coef_num_cell', 'coup_init1', '-v7.3');