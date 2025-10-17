CONFIG = SupplyNet_exp_config();

load(CONFIG.SAVEPATH_OUTPUTS);

plot_indices = 1:length(output.stats);

figure('Position', [0, 0, 500, 300]);

ha = tight_subplot(1, 1, [0, 0], [0.115, 0.025], [0.075, 0.020]);
axes(ha(1));
hold on;
box on;
grid on;

handle_LB = scatter([output.stats(plot_indices).support_size], [output.stats(plot_indices).DRO_LB], 25, 'blue', 'filled', ...
    'MarkerEdgeAlpha', 0.25, 'MarkerFaceAlpha', 0.25);
handle_UB = scatter([output.stats(plot_indices).support_size], [output.stats(plot_indices).DRO_UB], 25, 'red', 'filled', ...
    'MarkerEdgeAlpha', 0.25, 'MarkerFaceAlpha', 0.25);
set(gca, 'XLim', [0, 90]);

legend([handle_LB, handle_UB], ...
    '$\widehat{C}_{\mathrm{DRO}}^{\mathrm{LB}(t)}$', ...
    '$\widehat{C}_{\mathrm{DRO}}^{\mathrm{LB}(t)} + \hat{\varepsilon}_{\mathrm{sub}}^{(t)}$', ...
    'Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 15);

xlabel('$\big|\mathrm{supp}(\hat{\tau}^{(t)})\big|$', 'Interpreter', 'latex');
ylabel('objective');

figure('Position', [500, 0, 500, 300]);
ha = tight_subplot(1, 1, [0, 0], [0.115, 0.025], [0.075, 0.020]);
axes(ha(1));
hold on;
box on;
grid on;

handle_err = scatter([output.stats(plot_indices).support_size], [output.stats(plot_indices).error_sub], 25, 'black', 'filled', ...
    'MarkerEdgeAlpha', 0.25, 'MarkerFaceAlpha', 0.25);
set(gca, 'YScale', 'log');
set(gca, 'XLim', [0, 90]);

legend(handle_err, ...
    '$\hat{\varepsilon}_{\mathrm{sub}}^{(t)}$', ...
    'Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 15);

xlabel('$\big|\mathrm{supp}(\hat{\tau}^{(t)})\big|$', 'Interpreter', 'latex');
ylabel('sub-optimality');