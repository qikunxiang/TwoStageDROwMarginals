CONFIG = Scheduling_exp_config();

load(CONFIG.SAVEPATH_OUTPUTS);

figure('Position', [0, 0, 500, 300]);

ha = tight_subplot(1, 1, [0, 0], [0.100, 0.020], [0.070, 0.008]);
axes(ha(1));
hold on;
box on;

bar((1:length(st1_deci))', st1_deci, 'stacked', 'magenta', 'EdgeColor', 'none');

xlabel('task index');
ylabel('scheduled duration');