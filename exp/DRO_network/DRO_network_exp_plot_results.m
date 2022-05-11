load('exp/DRO_network/exp_rst_UB.mat', 'LSIP_UB_list');
load('exp/DRO_network/exp_rst_LB.mat', 'LB_errbar', ...
    'LB_samps_cell');


figure('Position', [100, 100, 250, 400]);
ha = tight_subplot(1, 1, [0.0, 0.0], [0.01, 0.010], [0.14, 0.025]);
hold on;

handle_UB = line([0, 2], [LSIP_UB_list(1), LSIP_UB_list(1)], ...
    'Color', 'blue', 'LineStyle', '--');

handle_LB = boxplot(LB_samps_cell{1}, 'Orientation', 'vertical', ...
    'Symbol', 'xr');

text(0.6, LSIP_UB_list(1) + 0.002, 'upper bound');
text(0.6, median(LB_samps_cell{1}), 'lower bound');

set(gca, 'YLim', [min(LB_samps_cell{1}) - 0.01, ...
    LSIP_UB_list(1) + 0.01]);
set(gca, 'XLim', [0.5, 1.2]);

set(gca, 'XTickLabel', ' ');

ylabel('worst-case expected cost');