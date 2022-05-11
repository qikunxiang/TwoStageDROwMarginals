load('exp/DRO_scheduling/exp_inputs.mat', 'knot_no_list', 'N', ...
    'err_bound_list', 'step_no');
load('exp/DRO_scheduling/exp_rst_UB.mat', 'LSIP_UB_list', ...
    'stage1_decision_cell');
load('exp/DRO_scheduling/exp_rst_LB.mat', ...
    'LB_errbar', 'LB_samps_cell');

xx = knot_no_list;

figure('Position', [100, 100, 400, 400]);
tight_subplot(1, 1, [0.0, 0.0], [0.08, 0.014], [0.080, 0.015]);

hold on;
LB_median = zeros(length(knot_no_list), 1);
LB_mean = zeros(length(knot_no_list), 1);

for stepid = 1:length(knot_no_list)
    LB_median(stepid) = median(LB_samps_cell{stepid});
    LB_mean(stepid) = mean(LB_samps_cell{stepid});
end

handle_UB = plot(xx, LSIP_UB_list, '-*', 'Color', 'blue');
set(gca, 'XScale', 'log');
xticklabel = get(gca, 'XTickLabel');
xtick = get(gca, 'XTick');
xtickmode = get(gca, 'XTickMode');
xticklabelmode = get(gca, 'XTickLabelMode');

boxplot(gca, horzcat(LB_samps_cell{:}), ...
    'Positions', xx, 'Widths', xx / xx(1) * 0.5, ...
    'Labels', xx, 'OutlierSize', 4, 'Symbol', 'xr');

handle_LB = plot(xx, LB_median, 'Marker', 'none', ...
    'LineStyle', '--', 'Color', 'red');

set(gca, 'XLim', [xx(1) * 0.92, xx(end) * 1.08]);
set(gca, 'XTick', xtick);
set(gca, 'XTickMode', xtickmode);
set(gca, 'XTickLabel', xticklabel);
set(gca, 'XTickLabelMode', xticklabelmode);
grid on;


legend([handle_UB, handle_LB], 'upper bound', 'lower bound');

xlabel('number of knots');
ylabel('worst-case expected total delay');



figure('Position', [500, 100, 400, 400]);
tight_subplot(1, 1, [0.0, 0.0], [0.08, 0.014], [0.10, 0.010]);
hold on;
bound_diff = LSIP_UB_list - LB_mean;

handle_theory = plot(xx, err_bound_list, '--^', 'Color', 'magenta');
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');

handle_computed = plot(xx, bound_diff, '-o', 'Color', 'black');

set(gca, 'XLim', [xx(1) * 0.92, xx(end) * 1.08]);
grid on;
box on;

legend([handle_theory, handle_computed], 'theoretical upper bound', ...
    'computed');

xlabel('number of knots');
ylabel('sub-optimality');


figure('Position', [900, 100, 400, 400]);
tight_subplot(1, 1, [0.0, 0.0], [0.08, 0.014], [0.085, 0.010]);
bar(1:N, stage1_decision_cell{step_no}, 'cyan');
xlabel('task index');
ylabel('scheduled execution time');
set(gca, 'XLim', [0.4, N + 0.6]);