CONFIG = MultiProdAssem_exp_config();

load(CONFIG.SAVEPATH_OUTPUTS);

figure('Position', [0, 0, 1000, 265]);

ha = tight_subplot(1, 4, [0, 0.024], [0.075, 0.095], [0.016, 0.004]);

Kpot = primal_sol.Kpot_duals_cell;

plot_index_list = [1; 2; 3; 4];

for i = 1:4
    axes(ha(i));

    plot_func = Kpot{plot_index_list(i)};
    
    line_x = repelem(plot_func(:, 1), 2, 1);
    line_y = reshape(plot_func(:, 3:4)', [], 1);
    
    plot(line_x, line_y, 'LineWidth', 1, 'Color', 'red');

    set(gca, 'XLim', [min(line_x) - 5, max(line_x) + 5]);
    set(gca, 'YLim', [min(line_y) - 0.2, max(line_y) + 0.2]);

    title(sprintf('$\\partial \\tilde{f}_{%d}$', plot_index_list(i)), 'Interpreter', 'latex', 'FontSize', 16);
    grid on;
end