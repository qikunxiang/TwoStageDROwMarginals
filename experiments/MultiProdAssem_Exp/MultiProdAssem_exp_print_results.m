CONFIG = MultiProdAssem_exp_config();

load(CONFIG.SAVEPATH_INPUTS);
load(CONFIG.SAVEPATH_OUTPUTS);

marginals = cell(product_num, 1);
marg_mean = zeros(product_num, 1);

for marg_id = 1:product_num
    marginals{marg_id} = eval(sprintf('%s(marginals_info{marg_id, 2}{:})', marginals_info{marg_id, 1}));
    marg_mean(marg_id) = marginals{marg_id}.computeMean();
end

products_degree = sum(assembly_parts > 0, 2);
parts_degree = sum(assembly_parts > 0, 1)';

fprintf('minimum mean demand = %.4f\n', min(marg_mean));
fprintf('maximum mean demand = %.4f\n', max(marg_mean));
fprintf('average mean demand = %.4f\n', mean(marg_mean));
fprintf('average number of parts needed for products = %.2f\n', mean(products_degree));
fprintf('average number of products used by parts = %.2f\n', mean(parts_degree));
fprintf('\n');

fprintf('error_tol = %.4e\n', tolerance);
fprintf('error_OT_tol = %.4e\n', OT_tolerance);
fprintf('\n');

fprintf('number of iterations = %4d\n', output.iter);
fprintf('DRO_LB = %.4f\n', DRO_LB);
fprintf('DRO_UB = %.4f\n', DRO_UB);
fprintf('error_sub = %.4e\n', error_sub);
fprintf('error_prob = %.4e\n', error_prob);
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