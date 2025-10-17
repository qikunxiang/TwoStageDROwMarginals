CONFIG = SupplyNet_exp_config();

load(CONFIG.SAVEPATH_INPUTS);
load(CONFIG.SAVEPATH_OUTPUTS);

G = digraph([network.S2P.sources; network.P2C.sources + supplier_num], ...
    [network.S2P.destinations + supplier_num; network.P2C.destinations + supplier_num + process_num]);

node_color = [repmat([1, 0, 0], supplier_num, 1); repmat([0, 0.75, 0], process_num, 1); repmat([0, 0, 1], customer_num, 1)];
line_width = [network.S2P.capacities; network.P2C.capacities];
edge_color = [2/3, 2/3, 2/3] .* (1 - [network.S2P.susceptibilities; network.P2C.susceptibilities]) ...
    + [1, 0, 0] .* [network.S2P.susceptibilities; network.P2C.susceptibilities];


figure('Position', [100, 100, 1000, 400]);
ha = tight_subplot(1, 1, [0.0, 0.0], [0.04, 0.004], [0.002, 0.002]);
hold on;
box on;

P = plot(G, '-', 'ArrowSize', 0, 'NodeLabel', {}, 'NodeColor', node_color, 'Marker', 'o', 'MarkerSize', 18, ...
    'LineWidth', line_width, 'EdgeColor', edge_color, 'EdgeAlpha', 0.5);


text(-4, 3, 'suppliers', 'FontSize', 20);
text(-4, 2.08, 'processing', 'FontSize', 20);
text(-4, 1.92, 'facilities', 'FontSize', 20);
text(-4, 1, 'customers', 'FontSize', 20);
text(-4, 0.24, 'first-stage', 'FontSize', 20);
text(-4, 0.08, 'decisions', 'FontSize', 20);

[process_node_x, process_node_sortind] = sort(P.XData(supplier_num + 1:supplier_num + process_num), 'ascend');
process_node_y = P.YData(supplier_num + 1:supplier_num + process_num);
process_node_y = process_node_y(process_node_sortind);

st1_deci = primal_sol.st1_deci(1:process_num);
st1_deci = st1_deci(process_node_sortind);

for proc_id = 1:process_num
    text(process_node_x(proc_id) - 0.26, process_node_y(proc_id), sprintf('%2d ', proc_id), 'FontSize', 12);
end

bar(process_node_x', st1_deci ./ process_capabilities_max * 0.75, 'FaceColor', [0, 0.75, 0]);


ax1 = gca;
ax1.XLim = [-4.4, 34.6];
ax1.YLim = [0, 3.10];
ax1.XTick = process_node_x;
ax1.XTickLabel = (1:process_num)';
ax1.YTick = [];
ax1.TickDir = 'none';