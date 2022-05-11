rng(1000, 'combRecursive');

N = 20;
K1 = 50;
K2ast = K1 + N;

% problem specification

% the mu and sigma^2 parameters of the mixture of truncated normal
% marginals 
mixtrnorm_mu_cell = cell(N, 1);
mixtrnorm_sig2_cell = cell(N, 1);
mixtrnorm_w_cell = cell(N, 1);
% truncation limits 
mixtrnorm_trunc = [0; 10];
mixtrnorm_compno = 3;

for i = 1:N
    mixtrnorm_mu_cell{i} = abs(randn(mixtrnorm_compno, 1) * 3);
    mixtrnorm_sig2_cell{i} = 1 ./ gamrnd(3, 1, mixtrnorm_compno, 1);
    mixtrnorm_w_cell{i} = ones(mixtrnorm_compno, 1) / mixtrnorm_compno;
end


while true
    assembly_mat = double(rand(N, K1) < 0.2);
    
    if all(sum(assembly_mat, 1) > 0) && all(sum(assembly_mat, 2) > 0)
        break;
    end
end

assembly_mat(assembly_mat == 1) = gamrnd(1, 1, nnz(assembly_mat), 1);
salvage_vec = gamrnd(0.5, 1, K1, 1);
price_vec = salvage_vec + gamrnd(1, 2, K1, 1);
production_cost_vec = assembly_mat * price_vec;
profit_vec = gamrnd(5, 5, N, 1);
revenue_vec = production_cost_vec + profit_vec;
stage2_ub = 1000 * ones(K2ast, 1);

DRO = DRO_multi_product_assembly(price_vec, assembly_mat, revenue_vec, ...
    salvage_vec, stage2_ub);

% approximation scheme
% the tolerance value used when solving LSIP
LSIP_tolerance = 1e-3;

% approximation scheme
knot_no_list = round(exp(linspace(log(5), log(100), 8)))';
step_no = length(knot_no_list);

knots_cell = cell(step_no, 1);
expval_cell = cell(step_no, 1);
coef_num_cell = cell(step_no, 1);

knots_hist_cell = cell(N, 1);

for i = 1:N
    % first generate the maximum number of knots
    [~, knots_hist_cell{i}, ~] = mixtruncnorm_momentset_greedy( ...
        mixtrnorm_mu_cell{i}, mixtrnorm_sig2_cell{i}, ...
        mixtrnorm_w_cell{i}, mixtrnorm_trunc, knot_no_list(end));
end

for step_id = 1:step_no
    knots = cell(N, 1);
    expval = cell(N, 1);
    coef_num_list = zeros(N, 1);
    
    for i = 1:N
        % retrive a subset of knots
        knots{i} = sort(knots_hist_cell{i}(1:knot_no_list(step_id)), ...
            'ascend');
        
        % compute the corresponding integrals
        expval{i} = mixtruncnorm_momentset(mixtrnorm_mu_cell{i}, ...
        mixtrnorm_sig2_cell{i}, mixtrnorm_w_cell{i}, knots{i});

        % the number of coefficients (with the constant intercept)
        coef_num_list(i) = length(knots{i});
    end
    
    knots_cell{step_id} = knots;
    expval_cell{step_id} = expval;
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

save('exp/DRO_assembly/exp_inputs.mat', ...
    'N', 'K1', 'K2ast', 'DRO', ...
    'mixtrnorm_mu_cell', 'mixtrnorm_sig2_cell', 'mixtrnorm_w_cell', ...
    'mixtrnorm_trunc', 'LSIP_tolerance', ...
    'knot_no_list', 'step_no', 'knots_cell', 'expval_cell', ...
    'coef_num_cell', 'coup_init1', '-v7.3');