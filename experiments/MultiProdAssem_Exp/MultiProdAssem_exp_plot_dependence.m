CONFIG = MultiProdAssem_exp_config();

load(CONFIG.SAVEPATH_OUTPUTS);

figure('Position', [0, 0, 800, 800]);

ha = tight_subplot(4, 4, [0.006, 0.006], [0.004, 0.027], [0.027, 0.004]);

plot_index_list = [1; 2; 3; 4];
scatter_list = (1:1e5)';

samp_min = min(wc_samples, [], 1)';
samp_max = max(wc_samples, [], 1)';

for i = 1:4
    for j = 1:4
        axes(ha((i - 1) * 4 + j));

        if i == j
            histogram(wc_samples(:, plot_index_list(i)), 80);
        else
            scatter(wc_samples(scatter_list, plot_index_list(j)), wc_samples(scatter_list, plot_index_list(i)), ...
                3, 'black', 'filled', 'MarkerEdgeAlpha', 0.1, 'MarkerFaceAlpha', 0.1);

            set(gca, 'YLim', [samp_min(plot_index_list(i)) - 0.1, samp_max(plot_index_list(i)) + 0.1]);
        end
    
        set(gca, 'XLim', [samp_min(plot_index_list(j)) - 0.1, samp_max(plot_index_list(j)) + 0.1]);
        set(gca, 'XTickLabel', []);
        set(gca, 'YTickLabel', []);
        set(gca, 'XAxisLocation', 'top');

        box on;

        if j == 1
            ylabel(sprintf('$\\mu_{%d}$', plot_index_list(i)), 'Interpreter', 'latex', 'FontSize', 18);
        end

        if i == 1
            xlabel(sprintf('$\\mu_{%d}$', plot_index_list(j)), 'Interpreter', 'latex', 'FontSize', 18);
        end
    end
end