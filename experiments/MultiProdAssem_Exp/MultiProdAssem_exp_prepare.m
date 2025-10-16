CONFIG = MultiProdAssem_exp_config();

rng(1200, 'combRecursive');

product_num = 50;
part_num = 100;


% problem specification
marginals_info = cell(product_num, 2);
marginals = cell(product_num, 1);

marginal_support_interval_lb = 2;
marginal_support_interval_ub = 10;

marginal_comp_num = 3;

marginal_options = struct;
marginal_options.bisect_bin_num = 2^10;
marginal_options.bisect_iter_num = 16;
marginal_options.boundary_prob_tolerance = 1e-15;

for marg_id = 1:product_num
    marginal_mean_list = abs(randn(marginal_comp_num, 1) * 3);
    marginal_std_list = sqrt(1 ./ gamrnd(3, 1, marginal_comp_num, 1));
    marginal_weight_list = ones(marginal_comp_num, 1) / marginal_comp_num;

    marginals_info{marg_id, 1} = 'TSDROMProbTruncMixNorm';
    marginals_info{marg_id, 2} = {marginal_support_interval_lb, marginal_support_interval_ub, ...
        marginal_weight_list, marginal_mean_list, marginal_std_list, marginal_options};
    marginals{marg_id} = eval(sprintf('%s(marginals_info{marg_id, 2}{:})', marginals_info{marg_id, 1}));
end

uncertain = TSDROMUncertain(marginals);

while true
    assembly_parts = double(rand(product_num, part_num) < 0.2);
    
    if all(sum(assembly_parts, 1) > 0) && all(sum(assembly_parts, 2) > 0)
        break;
    end
end

assembly_parts(assembly_parts == 1) = gamrnd(1, 1, nnz(assembly_parts), 1);
part_salvage_prices = gamrnd(0.5, 1, part_num, 1) + 0.1;
part_costs = part_salvage_prices + gamrnd(1, 2, part_num, 1);
production_costs = assembly_parts * part_costs;
product_profits = gamrnd(5, 5, product_num, 1);
product_revenues = production_costs + product_profits;

problem = TSDROMMultiProdAssem(uncertain, part_costs, assembly_parts, product_revenues, part_salvage_prices);

init_dual_point = problem.uncertain.computeOneDualPoint();

% An explicit choice of the initial dual point
% init_dual_point = [zeros(product_num, 1); assembly_parts * part_salvage_prices - product_revenues; ...
%     zeros(part_num, 1); -part_salvage_prices];

save(CONFIG.SAVEPATH_INPUTS, ...
    'part_num', ...
    'product_num', ...
    'marginals_info', ...
    'assembly_parts', ...
    'part_costs', ...
    'product_revenues', ...
    'part_salvage_prices', ...
    'init_dual_point', ...
    '-v7.3');
