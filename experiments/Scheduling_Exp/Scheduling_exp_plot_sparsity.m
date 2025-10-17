CONFIG = Scheduling_exp_config();

load(CONFIG.SAVEPATH_OUTPUTS);

figure('Position', [0, 0, 500, 300]);

ha = tight_subplot(1, 1, [0, 0], [0.115, 0.025], [0.077, 0.008]);
axes(ha(1));
hold on;
box on;
grid on;

handle_LB = scatter([output.stats.support_size], [output.stats.DRO_LB], 25, 'blue', 'filled', ...
    'MarkerEdgeAlpha', 0.4, 'MarkerFaceAlpha', 0.4);
handle_UB = scatter([output.stats.support_size], [output.stats.DRO_UB], 25, 'red', 'filled', ...
    'MarkerEdgeAlpha', 0.4, 'MarkerFaceAlpha', 0.4);
set(gca, 'XLim', [0, 130]);
set(gca, 'YLim', [-10000, 8000]);

legend([handle_LB, handle_UB], ...
    '$\widehat{C}_{\mathrm{DRO}}^{\mathrm{LB}(t)}$', ...
    '$\widehat{C}_{\mathrm{DRO}}^{\mathrm{LB}(t)} + \hat{\varepsilon}_{\mathrm{sub}}^{(t)}$', ...
    'Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 15);

xlabel('$\big|\mathrm{supp}(\hat{\tau}^{(t)})\big|$', 'Interpreter', 'latex');
ylabel('objective');

figure('Position', [500, 0, 500, 300]);
ha = tight_subplot(1, 1, [0, 0], [0.115, 0.025], [0.077, 0.008]);
axes(ha(1));
hold on;
box on;
grid on;

handle_err = scatter([output.stats.support_size], [output.stats.error_sub], 25, 'black', 'filled', ...
    'MarkerEdgeAlpha', 0.4, 'MarkerFaceAlpha', 0.4);
set(gca, 'YScale', 'log');
set(gca, 'XLim', [0, 130]);
set(gca, 'YLim', [5e-4, 5e4]);

legend(handle_err, ...
    '$\hat{\varepsilon}_{\mathrm{sub}}^{(t)}$', ...
    'Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 15);

xlabel('$\big|\mathrm{supp}(\hat{\tau}^{(t)})\big|$', 'Interpreter', 'latex');
ylabel('sub-optimality');