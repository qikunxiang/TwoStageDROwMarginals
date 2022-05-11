load('exp/DRO_assembly/exp_inputs.mat', 'knot_no_list');
load('exp/DRO_assembly/exp_rst_UB.mat', 'LSIP_UB_list');
load('exp/DRO_assembly/exp_rst_LB.mat', 'LB_errbar', ...
    'LB_samps_cell');

xx = knot_no_list;

LB_median = zeros(length(knot_no_list), 1);
LB_mean = zeros(length(knot_no_list), 1);

for stepid = 1:length(knot_no_list)
    LB_median(stepid) = median(LB_samps_cell{stepid});
    LB_mean(stepid) = mean(LB_samps_cell{stepid});
end

figure('Position', [100, 100, 800, 400]);
ha = tight_subplot(1, 2, [0.0, 0.07], [0.08, 0.010], [0.06, 0.005]);

axes(ha(1));
hold on;
handle_UB = plot(xx, LSIP_UB_list, '-*', 'Color', 'blue');
set(gca, 'XScale', 'log');

handle_LB = plot(xx, LB_median, 'Marker', 'x', ...
    'LineStyle', '--', 'Color', 'red');

set(gca, 'XLim', [xx(1) * 0.92, xx(end) * 1.08]);

box on;
grid on;

legend([handle_UB, handle_LB], 'upper bound', 'lower bound', ...
    'Location', 'southeast');

xlabel('number of knots');
ylabel('worst-case expected cost');

axes(ha(2));
hold on;

handle_UB = plot(xx(5:end), LSIP_UB_list(5:end), '-*', 'Color', 'blue');
set(gca, 'XScale', 'log');
xticklabel = get(gca, 'XTickLabel');
xtick = get(gca, 'XTick');
xtickmode = get(gca, 'XTickMode');
xticklabelmode = get(gca, 'XTickLabelMode');

boxplot(gca, horzcat(LB_samps_cell{5:end}), ...
    'Positions', xx(5:end), 'Widths', xx(5:end) / xx(5) * 1, ...
    'Labels', xx(5:end), 'OutlierSize', 4, 'Symbol', 'xr');

handle_LB = plot(xx(5:end), LB_median(5:end), 'Marker', 'none', ...
    'LineStyle', '--', 'Color', 'red');

set(gca, 'XLim', [xx(5) * 0.92, xx(end) * 1.08]);
set(gca, 'YLim', [min(LB_samps_cell{5}) - 0.5, LSIP_UB_list(5) + 0.5])
set(gca, 'XTick', xtick);
set(gca, 'XTickMode', xtickmode);
set(gca, 'XTickLabel', xticklabel);
set(gca, 'XTickLabelMode', xticklabelmode);

grid on;

legend([handle_UB, handle_LB], 'upper bound', 'lower bound', ...
    'Location', 'southeast');

xlabel('number of knots');
ylabel('worst-case expected cost');