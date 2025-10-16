CONFIG = Scheduling_exp_config();

rng(7000, 'combRecursive');

% problem specification
task_num = 100;
time_max = 200;

delay_weights = ones(task_num, 1) + rand(task_num, 1);
delay_weights(end) = delay_weights(end) * 10;

marginal_support_interval_lb = 0;
marginal_support_interval_ub = 5;
marginal_weight_list = [0.4; 0.4; 0.2];
marginal_mean_list = [0.3; 1; 2];
marginal_std_list = [0.2; 0.4; 1];

marginal_options = struct;
marginal_options.bisect_bin_num = 2^10;
marginal_options.bisect_iter_num = 16;
marginal_options.boundary_prob_tolerance = 1e-15;

marginals_info = cell(task_num, 2);
marginals = cell(task_num, 1);

for marg_id = 1:task_num
    marginals_info{marg_id, 1} = 'TSDROMProbTruncMixNorm';
    marginals_info{marg_id, 2} = {marginal_support_interval_lb, marginal_support_interval_ub, ...
        marginal_weight_list, marginal_mean_list, marginal_std_list, marginal_options};
    marginals{marg_id} = eval(sprintf('%s(marginals_info{marg_id, 2}{:})', marginals_info{marg_id, 1}));
end

uncertain = TSDROMUncertain(marginals);

problem = TSDROMScheduling(uncertain, time_max, delay_weights);

init_dual_point = problem.uncertain.computeOneDualPoint();

save(CONFIG.SAVEPATH_INPUTS, ...
    'task_num', ...
    'time_max', ...
    'marginals_info', ...
    'delay_weights', ...
    'init_dual_point', ...
    '-v7.3');
