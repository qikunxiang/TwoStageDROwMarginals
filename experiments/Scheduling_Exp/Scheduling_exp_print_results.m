CONFIG = Scheduling_exp_config();

load(CONFIG.SAVEPATH_INPUTS);
load(CONFIG.SAVEPATH_OUTPUTS);

marginals = cell(task_num, 1);

for marg_id = 1:task_num
    marginals{marg_id} = eval(sprintf('%s(marginals_info{marg_id, 2}{:})', marginals_info{marg_id, 1}));
end

fprintf('mean duration = %.4f\n', marginals{1}.computeMean());
fprintf('\n');

fprintf('error_tol = %.4e\n', tolerance);
fprintf('error_OT_tol = %.4e\n', OT_tolerance);
fprintf('\n');

fprintf('number of iterations = %4d\n', output.iter);
fprintf('DRO_LB = %.4f\n', DRO_LB);
fprintf('DRO_UB = %.4f\n', DRO_UB);
fprintf('error_sub = %.4e\n', error_sub);
fprintf('error_prob = %.4e\n', error_prob);
fprintf('error_sub_rel = %.4e\n', error_sub / abs(DRO_UB));
fprintf('error_prob_rel = %.4e\n', error_prob / abs(DRO_UB));
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

for task_id = 1:task_num
    [~, uind, umap] = unique(round(dual_sol.dual_points(:, task_id), 4));

    uprob = accumarray(umap, dual_sol.dual_points_probs);

    [uprob, sorted_ind] = sort(uprob, 'descend');
    uind = uind(sorted_ind);
    
    fprintf('Marginal %d, %2d atoms:\n', task_id, length(uind));

    for a_id = 1:length(uind)
        fprintf('%15.4f  ', dual_sol.dual_points(uind(a_id), task_id));
    end
    
    fprintf('\n');

    for a_id = 1:length(uind)
        fprintf('%15.4e  ', uprob(a_id));
    end

    fprintf('\n');

end