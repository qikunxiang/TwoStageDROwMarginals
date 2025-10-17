CONFIG = SupplyNet_exp_config();

load(CONFIG.SAVEPATH_INPUTS);
load(CONFIG.SAVEPATH_OUTPUTS);

demand_marginals = cell(customer_num, 1);
demand_mean = zeros(customer_num, 1);

for marg_id = 1:customer_num
    demand_marginals{marg_id} = eval(sprintf('%s(demand_marginals_info{marg_id, 2}{:})', demand_marginals_info{marg_id, 1}));
    demand_mean(marg_id) = demand_marginals{marg_id}.computeMean();
end

failure_probs = zeros(edge_susceptible_num, 1);

for edge_id = 1:edge_susceptible_num
    failure_probs(edge_id) = edge_failure_marginals_info{edge_id, 2}{1};
end

fprintf('minimum mean demand = %.4f\n', min(demand_mean));
fprintf('maximum mean demand = %.4f\n', max(demand_mean));
fprintf('average mean demand = %.4f\n', mean(demand_mean));
fprintf('minimum failure probability = %.4f\n', min(failure_probs));
fprintf('maximum failure probability = %.4f\n', max(failure_probs));
fprintf('average failure probability = %.4f\n', mean(failure_probs));
fprintf('\n');

fprintf('error_tol = %.4e\n', tolerance);
fprintf('error_OT_tol = %.4e\n', OT_tolerance);
fprintf('\n');

fprintf('number of iterations = %4d\n', output.iter);
fprintf('DRO_LB = %.4f\n', DRO_LB);
fprintf('DRO_UB = %.4f\n', DRO_UB);
fprintf('error_sub = %.4f\n', error_sub);
fprintf('error_prob = %.4f\n', error_prob);
fprintf('error_sub_rel = %.4f%%\n', error_sub / abs(DRO_UB) * 100);
fprintf('error_prob_rel = %.4f%%\n', error_prob / abs(DRO_UB) * 100);
fprintf('\n');

Kpot_knot_num = zeros(length(primal_sol.Kpot_duals_cell), 1);

for marg_id = 1:length(primal_sol.Kpot_duals_cell)
    Kpot_knot_num(marg_id) = size(primal_sol.Kpot_duals_cell{marg_id}, 1);
end

fprintf('minimum knot number = %4d\n', min(Kpot_knot_num));
fprintf('maximum knot number = %4d\n', max(Kpot_knot_num));
fprintf('average knot number = %4.2f\n', mean(Kpot_knot_num));
fprintf('maximum support size = %4d\n', max([output.stats.support_size]));
fprintf('support size = %4d\n', size(dual_sol.dual_points, 1));
fprintf('\n');