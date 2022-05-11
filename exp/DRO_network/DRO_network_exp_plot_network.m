load('exp/DRO_network/exp_inputs.mat', 'edges', 'supplier_no', ...
    'processing_no', 'customer_no', 'supply_vec', 'failure_prob',...
    'max_process_capa_vec', 'process_invest_cost_vec', ...
    'knots_cell', 'expval_cell');
load('exp/DRO_network/exp_rst_UB.mat', 'stage1_decision_cell');

edge_s2p_no = length(edges.s2p.from);
edge_p2c_no = length(edges.p2c.from);

figure('Position', [100, 100, 800, 400]);
ha = tight_subplot(1, 1, [0.0, 0.0], [0.01, 0.01], [0.01, 0.01]);

hold on;

G = digraph([edges.s2p.from; edges.p2c.from + supplier_no], ...
    [edges.s2p.to + supplier_no; edges.p2c.to ...
    + supplier_no + processing_no]);

linestyle = cell(edge_s2p_no + edge_p2c_no, 1);

for i = 1:edge_s2p_no
    linestyle{i} = '-';
end

for i = 1:edge_p2c_no
    linestyle{edge_s2p_no + i} = '-';
end

edgewidth = [edges.s2p.capacity; edges.p2c.capacity];

edgealpha = [edges.s2p.cost; edges.p2c.cost];
edgealpha = (edgealpha / max(edgealpha)) .^ 0.35;

edgecolor = zeros(edge_s2p_no + edge_p2c_no, 3);
edgecolor([edges.s2p.susceptible; edges.p2c.susceptible], :) ...
    = repmat([1, 0, 0], sum([edges.s2p.susceptible; ...
    edges.p2c.susceptible]), 1);

edgecolor = ones(edge_s2p_no + edge_p2c_no, 3) .* (1 - edgealpha) ...
    + edgecolor .* edgealpha;

marker = cell(supplier_no + processing_no + customer_no, 1);

for i = 1:supplier_no
    marker{i} = 'o';
end

for i = 1:processing_no
    marker{supplier_no + i} = 'o';
end

for i = 1:customer_no
    marker{supplier_no + processing_no + i} = 'o';
end

mean_demand = zeros(customer_no, 1);

for i = 1:customer_no
    mean_demand(i) = knots_cell{1}{i}' * expval_cell{1}{i};
end

markersize_s = supply_vec / max(supply_vec) * 20;
markersize_p = max_process_capa_vec / max(max_process_capa_vec) * 20;
markersize_p2 = markersize_p / 20 * 16;
markersize_p3 = stage1_decision_cell{1}(1:processing_no) ...
    / max(max(max_process_capa_vec)) * 16;
markersize_c = mean_demand / max(mean_demand) * 20;
markersize = [markersize_s; markersize_p; markersize_c];
markersize2 = [markersize_s; markersize_p2; markersize_c];
markersize3 = [markersize_s; markersize_p3; markersize_c];

process_invest_cost_normalized = process_invest_cost_vec ...
    / max(process_invest_cost_vec);

nodecolor = [repmat([1, 0, 1], supplier_no, 1);
    repmat([0, 0, 1], processing_no, 1) .* process_invest_cost_normalized ...
    + repmat([1, 1, 1], processing_no, 1) ...
    .* (1 - process_invest_cost_normalized);
    repmat([0, 1, 0], customer_no, 1)];

nodecolor2 = [repmat([1, 0, 0], supplier_no, 1);
    repmat([1, 1, 1], processing_no, 1);
    repmat([0, 1, 0], customer_no, 1)];

nodecolor3 = [repmat([1, 0, 1], supplier_no, 1);
    repmat([0, 0, 1], processing_no, 1);
    repmat([0, 1, 0], customer_no, 1)];

edgelabel = cell(edge_s2p_no + edge_p2c_no, 1);

susc_counter = 0;
for j = 1:edge_s2p_no
    if edges.s2p.susceptible(j)
        susc_counter = susc_counter + 1;
        edgelabel{j} = sprintf('%d%%', round(failure_prob(susc_counter) * 100));

        if susc_counter == 1
            edgelabel{j} = [edgelabel{j}, '                          '];
        elseif susc_counter == 2
            edgelabel{j} = [edgelabel{j}, '                          '];
        elseif susc_counter == 6
            edgelabel{j} = [edgelabel{j}, '                          '];
        elseif susc_counter == 9
            edgelabel{j} = ['                          ', edgelabel{j}];
        elseif susc_counter == 14
            edgelabel{j} = ['                                    ', ...
                edgelabel{j}];
        elseif susc_counter == 10
            edgelabel{j} = ['                                    ', ...
                edgelabel{j}];
        end
    else
        edgelabel{j} = '';
    end
end

for j = 1:edge_p2c_no
    if edges.p2c.susceptible(j)
        susc_counter = susc_counter + 1;
        edgelabel{edge_s2p_no + j}...
            = sprintf('%d%%', round(failure_prob(susc_counter) * 100));

        if susc_counter == 23
            edgelabel{edge_s2p_no+ j} = ['                          ', ...
                edgelabel{edge_s2p_no + j}];
        elseif susc_counter == 24
            edgelabel{edge_s2p_no + j} = [edgelabel{edge_s2p_no + j}, ...
                '                          '];
        end
    else
        edgelabel{edge_s2p_no + j} = '';
    end
end


plot(G, 'LineStyle', linestyle, 'ArrowSize', 0, 'NodeLabel', {}, ...
    'NodeColor', nodecolor, 'Marker', marker, ...
    'MarkerSize', markersize, 'LineWidth', edgewidth, ...
    'EdgeColor', edgecolor, 'EdgeAlpha', 1);

plot(G, 'LineStyle', 'none', 'ArrowSize', 0, 'NodeLabel', {}, ...
    'NodeColor', nodecolor2, 'Marker', marker, ...
    'MarkerSize', markersize2, 'LineWidth', edgewidth, ...
    'EdgeColor', edgecolor, 'EdgeAlpha', 1);

plot(G, 'LineStyle', 'none', 'ArrowSize', 0, 'NodeLabel', {}, ...
    'NodeColor', nodecolor3, 'Marker', marker, ...
    'MarkerSize', markersize3, 'LineWidth', edgewidth, ...
    'EdgeColor', edgecolor, 'EdgeAlpha', 1, 'EdgeLabel', edgelabel);

text(-4, 3, 'suppliers', 'FontSize', 20);
text(-4, 2.08, 'processing', 'FontSize', 20);
text(-4, 1.92, 'facilities', 'FontSize', 20);
text(-4, 1, 'customers', 'FontSize', 20);

set(gca, 'XLim', [-4.5, 26]);
set(gca, 'YLim', [0.9, 3.1]);
ax1 = gca;   
ax1.YAxis.Visible = 'off';
ax1.XAxis.Visible = 'off';