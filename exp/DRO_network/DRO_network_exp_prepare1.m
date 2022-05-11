rng(1000, 'combRecursive');

supplier_no = 15;
processing_no = 20;
customer_no = 10;

edge_s2p_no = 90;
edge_p2c_no = 60;
edge_s2p_scp_no = 15;
edge_p2c_scp_no = 10;

while true

    edges = struct;
    edges.s2p = struct;
    edges.s2p.from = repelem((1:supplier_no)', processing_no, 1);
    edges.s2p.to = repmat((1:processing_no)', supplier_no, 1);
    
    % randomly keep some of the edges
    s2p_keep_list = randsample(supplier_no * processing_no, ...
        edge_s2p_no, false);
    edges.s2p.from = edges.s2p.from(s2p_keep_list);
    edges.s2p.to = edges.s2p.to(s2p_keep_list);

    % check if all suppliers and processing facilities are connected
    if length(unique(edges.s2p.from)) ~= supplier_no ...
            || length(unique(edges.s2p.to)) ~= processing_no
        % if not, re-generate the edges
        continue;
    end

    edges.s2p.cost = 1 ./ gamrnd(3, 0.5, edge_s2p_no, 1);
    edges.s2p.capacity = 1 ./ gamrnd(5, 0.4, edge_s2p_no, 1);

    % set the cheapest edges to susceptible
    [~, edges_s2p_order] = sort(edges.s2p.cost, 'ascend');
    edges.s2p.susceptible = false(edge_s2p_no, 1);
    edges.s2p.susceptible(edges_s2p_order(1:edge_s2p_scp_no)) = true;

    
    edges.p2c = struct;
    edges.p2c.from = repelem((1:processing_no)', customer_no, 1);
    edges.p2c.to = repmat((1:customer_no)', processing_no, 1);

    % randomly keep some of the edges
    p2c_keep_list = randsample(processing_no * customer_no, ...
        edge_p2c_no, false);
    edges.p2c.from = edges.p2c.from(p2c_keep_list);
    edges.p2c.to = edges.p2c.to(p2c_keep_list);

    % check if all processing facilities and customers are connected
    if length(unique(edges.p2c.from)) ~= processing_no ...
            || length(unique(edges.p2c.to)) ~= customer_no
        % if not, re-generate the edges
        continue;
    end

    edges.p2c.cost = 1 ./ gamrnd(3, 0.5, edge_p2c_no, 1);
    edges.p2c.capacity = 1 ./ gamrnd(5, 0.4, edge_p2c_no, 1);

    % set the cheapest edges to susceptible
    [~, edges_p2c_order] = sort(edges.p2c.cost, 'ascend');
    edges.p2c.susceptible = false(edge_p2c_no, 1);
    edges.p2c.susceptible(edges_p2c_order(1:edge_p2c_scp_no)) = true;
    
    edge_no = length(edges.s2p.from) + length(edges.p2c.from);
    edge_scp_no = sum(edges.s2p.susceptible) + sum(edges.p2c.susceptible);
    
    max_demand_vec = ones(customer_no, 1) * 2;
    supply_vec = rand(supplier_no, 1);
    supply_vec = supply_vec / sum(supply_vec) * sum(max_demand_vec) * 1.2;
    max_process_capa_vec = ones(processing_no, 1) / processing_no ...
        * sum(max_demand_vec) * 2;
    process_invest_cost_vec = 1 ./ gamrnd(3, 0.2, processing_no, 1);
    
    balance_dual_lb = -1000 * ones(processing_no, 1);
    balance_dual_ub = 1000 * ones(processing_no, 1);
    supply_dual_lb = -1000 * ones(supplier_no, 1);
    demand_dual_ub = 1000 * ones(customer_no, 1);
    processing_dual_lb = -1000 * ones(processing_no, 1);
    capacity_dual_lb = -1000 * ones(edge_no, 1);
    
    
    % problem specification
    
    % the mu and sigma^2 parameters of the mixture of truncated normal
    % marginals 
    mixtrnorm_mu_cell = cell(customer_no, 1);
    mixtrnorm_sig2_cell = cell(customer_no, 1);
    mixtrnorm_w_cell = cell(customer_no, 1);
    % truncation limits 
    mixtrnorm_trunc_cell = cell(customer_no, 1);
    mixtrnorm_compno = 3;
    
    for i = 1:customer_no
        mixtrnorm_trunc_cell{i} = [0; max_demand_vec(i)];
        mixtrnorm_mu_cell{i} = abs(randn(mixtrnorm_compno, 1) + 1);
        mixtrnorm_sig2_cell{i} = 1 ./ gamrnd(3, 1, mixtrnorm_compno, 1);
        mixtrnorm_w_cell{i} = ones(mixtrnorm_compno, 1) / mixtrnorm_compno;
    end
    
    failure_prob = betarnd(1, 9, edge_scp_no, 1);
    
    
    
    DRO = DRO_supply_chain_network(edges, supply_vec, ...
        max_process_capa_vec, process_invest_cost_vec, ...
        max_demand_vec, balance_dual_lb, balance_dual_ub, ...
        supply_dual_lb, demand_dual_ub, processing_dual_lb, ...
        capacity_dual_lb);
    SP = SP_supply_chain_network(edges, supply_vec, ...
        max_process_capa_vec, process_invest_cost_vec, max_demand_vec);

    try
        SP_approx(SP, [max_demand_vec; ones(edge_scp_no, 1)]', 1);
    catch exception
        warning('the problem is infeasible, retry...');
        continue;
    end

    break;
end

% approximation scheme
% the tolerance value used when solving LSIP
LSIP_tolerance = 1e-2;

% approximation scheme
knot_no_list = 40;
step_no = length(knot_no_list);

knots_cell = cell(step_no, 1);
expval_cell = cell(step_no, 1);
coef_num_cell = cell(step_no, 1);

knots_hist_cell = cell(customer_no, 1);

for i = 1:customer_no
    % first generate the maximum number of knots
    [~, knots_hist_cell{i}, ~] = mixtruncnorm_momentset_greedy( ...
        mixtrnorm_mu_cell{i}, mixtrnorm_sig2_cell{i}, ...
        mixtrnorm_w_cell{i}, mixtrnorm_trunc_cell{i}, knot_no_list(end));
end

N = customer_no + edge_scp_no;
K1 = processing_no + edge_no;
K2ast = processing_no + supplier_no + customer_no + processing_no ...
    + edge_no;

for step_id = 1:step_no
    knots = cell(N, 1);
    expval = cell(N, 1);
    coef_num_list = zeros(N, 1);
    
    for i = 1:customer_no
        % retrive a subset of knots
        knots{i} = sort(knots_hist_cell{i}(1:knot_no_list(step_id)), ...
            'ascend');
        
        % compute the corresponding integrals
        expval{i} = mixtruncnorm_momentset(mixtrnorm_mu_cell{i}, ...
        mixtrnorm_sig2_cell{i}, mixtrnorm_w_cell{i}, knots{i});

        % the number of coefficients (with the constant intercept)
        coef_num_list(i) = length(knots{i});
    end

    for i = 1:edge_scp_no
        knots{customer_no + i} = [0; 1];
        expval{customer_no + i} = [1 - failure_prob(i); failure_prob(i)];
        coef_num_list(customer_no + i) = 2;
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

save('exp/DRO_network/exp_inputs.mat', ...
    'N', 'K1', 'K2ast', 'DRO', ...
    'supplier_no', 'processing_no', 'customer_no', 'edges', ...
    'edge_no', 'edge_scp_no', 'supply_vec', 'max_process_capa_vec', ...
    'process_invest_cost_vec', 'max_demand_vec', 'failure_prob', ...
    'mixtrnorm_mu_cell', 'mixtrnorm_sig2_cell', 'mixtrnorm_w_cell', ...
    'mixtrnorm_trunc_cell', 'LSIP_tolerance', ...
    'knot_no_list', 'step_no', 'knots_cell', 'expval_cell', ...
    'coef_num_cell', 'coup_init1', '-v7.3');