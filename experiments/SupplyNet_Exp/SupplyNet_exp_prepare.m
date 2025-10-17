CONFIG = SupplyNet_exp_config();

rng(1000, 'combRecursive');

% problem specification
supplier_num = 20;
process_num = 20;
customer_num = 20;

edge_S2P_num = 100;
edge_P2C_num = 100;
edge_S2P_susceptible_num = 20;
edge_P2C_susceptible_num = 20;

while true
    network = struct;
    network.S2P = struct;
    network.S2P.sources = repelem((1:supplier_num)', process_num, 1);
    network.S2P.destinations = repmat((1:process_num)', supplier_num, 1);
    
    % randomly keep some of the edges
    S2P_keep_list = randsample(supplier_num * process_num, edge_S2P_num, false);
    network.S2P.sources = network.S2P.sources(S2P_keep_list);
    network.S2P.destinations = network.S2P.destinations(S2P_keep_list);

    % check if all suppliers and processing facilities are connected
    if length(unique(network.S2P.sources)) ~= supplier_num || length(unique(network.S2P.destinations)) ~= process_num
        % if not, re-generate the edges
        continue;
    end

    network.S2P.costs = 1 ./ gamrnd(3, 0.5, edge_S2P_num, 1);
    network.S2P.capacities = 1 ./ gamrnd(5, 0.2, edge_S2P_num, 1);

    % set the cheapest edges to susceptible
    [~, edges_S2P_order] = sort(network.S2P.costs, 'ascend');
    network.S2P.susceptibilities = false(edge_S2P_num, 1);
    network.S2P.susceptibilities(edges_S2P_order(1:edge_S2P_susceptible_num)) = true;

    
    network.P2C = struct;
    network.P2C.sources = repelem((1:process_num)', customer_num, 1);
    network.P2C.destinations = repmat((1:customer_num)', process_num, 1);

    % randomly keep some of the edges
    P2C_keep_list = randsample(process_num * customer_num, edge_P2C_num, false);
    network.P2C.sources = network.P2C.sources(P2C_keep_list);
    network.P2C.destinations = network.P2C.destinations(P2C_keep_list);

    % check if all processing facilities and customers are connected
    if length(unique(network.P2C.sources)) ~= process_num || length(unique(network.P2C.destinations)) ~= customer_num
        % if not, re-generate the edges
        continue;
    end

    network.P2C.costs = 1 ./ gamrnd(3, 0.5, edge_P2C_num, 1);
    network.P2C.capacities = 1 ./ gamrnd(5, 0.2, edge_P2C_num, 1);

    % set the cheapest edges to susceptible
    [~, edges_P2C_order] = sort(network.P2C.costs, 'ascend');
    network.P2C.susceptibilities = false(edge_P2C_num, 1);
    network.P2C.susceptibilities(edges_P2C_order(1:edge_P2C_susceptible_num)) = true;
    
    edge_num = edge_S2P_num + edge_P2C_num;
    edge_susceptible_num = edge_S2P_susceptible_num + edge_P2C_susceptible_num;
    
    demand_max = ones(customer_num, 1) * 2;
    supplies = rand(supplier_num, 1);
    supplies = supplies / sum(supplies) * sum(demand_max) * 1.2;
    process_capabilities_max = ones(process_num, 1) / process_num * sum(demand_max) * 2;
    process_investment_costs = 1 ./ gamrnd(3, 0.2, process_num, 1);
    
    % the mixture of truncated normal marginals 
    demand_marginals_info = cell(customer_num, 2);
    demand_marginals = cell(customer_num, 1);

    demand_marginal_comp_num = 3;

    demand_marginal_options = struct;
    demand_marginal_options.bisect_bin_num = 2^10;
    demand_marginal_options.bisect_iter_num = 16;
    demand_marginal_options.boundary_prob_tolerance = 1e-15;
    
    for cust_id = 1:customer_num
        demand_marginal_support_interval_lb = demand_max(cust_id) / 2;
        demand_marginal_support_interval_ub = demand_max(cust_id);
        demand_marginal_weight_list = ones(demand_marginal_comp_num, 1) / demand_marginal_comp_num;
        demand_marginal_mean_list = abs(randn(demand_marginal_comp_num, 1) + 1);
        demand_marginal_std_list = 1 ./ sqrt(gamrnd(3, 1, demand_marginal_comp_num, 1));

        demand_marginals_info{cust_id, 1} = 'TSDROMProbTruncMixNorm';
        demand_marginals_info{cust_id, 2} = {demand_marginal_support_interval_lb, demand_marginal_support_interval_ub, ...
            demand_marginal_weight_list, demand_marginal_mean_list, demand_marginal_std_list, demand_marginal_options};
        demand_marginals{cust_id} = eval(sprintf('%s(demand_marginals_info{cust_id, 2}{:})', demand_marginals_info{cust_id, 1}));
    end
    

    edge_failure_probabilities = betarnd(1, 9, edge_susceptible_num, 1);
    edge_failure_marginals_info = cell(edge_susceptible_num, 2);
    edge_failure_marginals = cell(edge_susceptible_num, 1);

    for edge_id = 1:edge_susceptible_num
        edge_failure_marginals_info{edge_id, 1} = 'TSDROMProbBernoulli';
        edge_failure_marginals_info{edge_id, 2} = {edge_failure_probabilities(edge_id)};
        edge_failure_marginals{edge_id} = eval(sprintf('%s(edge_failure_marginals_info{edge_id, 2}{:})', ...
            edge_failure_marginals_info{edge_id, 1}));
    end
    
    uncertain = TSDROMSupplyNetUncertain(demand_marginals, edge_failure_marginals);
    problem = TSDROMSupplyNet(uncertain, network, supplies, process_capabilities_max, process_investment_costs);
    
    if ~problem.checkFeasibility()
        warning('the problem is infeasible, retry...');
        continue;
    end

    break;
end

init_dual_point = problem.uncertain.computeOneDualPoint();

save(CONFIG.SAVEPATH_INPUTS, ...
    'supplier_num', ...
    'process_num', ...
    'customer_num', ...
    'edge_S2P_num', ...
    'edge_P2C_num', ...
    'edge_num', ...
    'edge_S2P_susceptible_num', ...
    'edge_P2C_susceptible_num', ...
    'edge_susceptible_num', ...
    'demand_marginals_info', ...
    'edge_failure_marginals_info', ...
    'network', ...
    'supplies', ...
    'demand_max', ...
    'process_capabilities_max', ...
    'process_investment_costs', ...
    'init_dual_point', ...
    '-v7.3');
