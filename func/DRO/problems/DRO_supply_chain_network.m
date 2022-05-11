function DRO = DRO_supply_chain_network(edges, supply_vec, ...
    max_process_capa_vec, process_invest_cost_vec, max_demand_vec, ...
    balance_dual_lb, balance_dual_ub, supply_dual_lb, demand_dual_ub, ...
    processing_dual_lb, capacity_dual_lb)
% Construct the DRO problem specification for a supply chain network design
% problem with edge failure
% Inputs: 
%       edges: struct containing information about the edges in the supply
%       chain network
%           s2p: struct containing information about edges from suppliers
%           to processing facilities
%               from: vector of supplier indices
%               to: vector of processing facility indices
%               cost: vector of transportation costs
%               capacity: vector of maximum capacities
%               susceptible: boolean vector indicating whether the edge is
%               susceptible to failure
%           p2c: struct containing information about edges from processing
%           facilities to customers; fields are same as above
%       supply_vec: vector containing the supply of each supplier (sum of
%       the supply must exceed the sum of the maximum demand)
%       max_process_capa_vec: vector containing the full processing
%       capabilities of the processing facilities
%       process_invest_cost_vec: the per unit cost of investment into the
%       processing capability of each processing facility
%       max_demand_vec: vector containing the maximum demand of each
%       customer
%       balance_dual_lb: lower bounds for the dual variables of the balance
%       constraints
%       balance_dual_ub: upper bounds for the dual variables of the balance
%       constraints
%       supply_dual_lb: lower bounds for the dual variables of the supply
%       constraints
%       demand_dual_ub: upper bounds for the dual variables of the demand
%       constraints
%       processing_dual_lb: lower bounds for the dual variables of the
%       processing capability constraints
%       capacity_dual_lb: lower bounds for the dual variables of the
%       capacity constraints
% Output:
%       DRO: struct containing the DRO specifications

supplier_no = length(supply_vec);
processing_no = length(max_process_capa_vec);
customer_no = length(max_demand_vec);
edge_s2p_no = length(edges.s2p.from);
edge_p2c_no = length(edges.p2c.from);
edge_no = edge_s2p_no + edge_p2c_no;
edge_s2p_scp_no = sum(edges.s2p.susceptible);
edge_p2c_scp_no = sum(edges.p2c.susceptible);
edge_scp_no = edge_s2p_scp_no + edge_p2c_scp_no;

% check the inputs
assert(sum(supply_vec) >= sum(max_demand_vec), 'supply is too low');
assert(sum(max_process_capa_vec) >= sum(max_demand_vec), ...
    'processing capability is too low');
assert(length(process_invest_cost_vec) == processing_no);

assert(length(edges.s2p.to) == edge_s2p_no ...
    && length(edges.s2p.cost) == edge_s2p_no ...
    && length(edges.s2p.capacity) == edge_s2p_no ...
    && length(edges.s2p.susceptible) == edge_s2p_no ...
    && length(edges.p2c.to) == edge_p2c_no ...
    && length(edges.p2c.cost) == edge_p2c_no ...
    && length(edges.p2c.capacity) == edge_p2c_no ...
    && length(edges.p2c.susceptible) == edge_p2c_no, ...
    'edge info invalid');

assert(length(unique(edges.s2p.from)) == supplier_no ...
    && all(unique(edges.s2p.from) == (1:supplier_no)') ...
    && length(unique(edges.s2p.to)) == processing_no ...
    && all(unique(edges.s2p.to) == (1:processing_no)') ...
    && length(unique(edges.p2c.from)) == processing_no ...
    && all(unique(edges.p2c.from) == (1:processing_no)') ...
    && length(unique(edges.p2c.to)) == customer_no ...
    && all(unique(edges.p2c.to) == (1:customer_no)'), ...
    'some notes are not connected');

assert(length(balance_dual_lb) == processing_no);
assert(length(balance_dual_ub) == processing_no);
assert(length(supply_dual_lb) == supplier_no);
assert(length(demand_dual_ub) == customer_no);
assert(length(processing_dual_lb) == processing_no);
assert(length(capacity_dual_lb) == edge_no);


from_s_cell = cell(supplier_no, 1);
to_p_cell = cell(processing_no, 1);
from_p_cell = cell(processing_no, 1);
to_c_cell = cell(customer_no, 1);

for s = 1:supplier_no
    from_s_cell{s} = double(sparse(edges.s2p.from == s))';
end

for p = 1:processing_no
    to_p_cell{p} = double(sparse(edges.s2p.to == p))';
    from_p_cell{p} = double(sparse(edges.p2c.from == p))';
end

for c = 1:customer_no
    to_c_cell{c} = double(sparse(edges.p2c.to == c))';
end

from_s_mat = vertcat(from_s_cell{:});
to_p_mat = vertcat(to_p_cell{:});
from_p_mat = vertcat(from_p_cell{:});
to_c_mat = vertcat(to_c_cell{:});

p_in_out_mat = [to_p_mat, -from_p_mat];
s_out_mat = [from_s_mat, sparse(supplier_no, edge_p2c_no)];
c_in_mat = [sparse(customer_no, edge_s2p_no), to_c_mat];
p_in_mat = [to_p_mat, sparse(processing_no, edge_p2c_no)];

s2p_scp_list = edges.s2p.susceptible == 1;
p2c_scp_list = edges.p2c.susceptible == 1;
capacity_const_vec = [edges.s2p.capacity; edges.p2c.capacity];

capacity_link_mat_s2p = sparse(find(s2p_scp_list), 1:edge_s2p_scp_no, ...
    edges.s2p.capacity(s2p_scp_list), edge_s2p_no, edge_s2p_scp_no);
capacity_link_mat_p2c = sparse(find(p2c_scp_list), 1:edge_p2c_scp_no, ...
    edges.p2c.capacity(p2c_scp_list), edge_p2c_no, edge_p2c_scp_no);
capacity_link_mat = blkdiag(capacity_link_mat_s2p, capacity_link_mat_p2c);

K2ast = processing_no * 2 + supplier_no + customer_no + edge_no;

DRO = struct;
DRO.stage1 = struct;
DRO.stage1.c1 = [process_invest_cost_vec; zeros(edge_no, 1)];
DRO.stage1.Lin = [[sparse(supplier_no + customer_no, processing_no); ...
    -speye(processing_no)], [s_out_mat; -c_in_mat; p_in_mat]];
DRO.stage1.qin = [supply_vec; -max_demand_vec; zeros(processing_no, 1)];
DRO.stage1.Leq = [sparse(processing_no, processing_no), p_in_out_mat];
DRO.stage1.qeq = zeros(processing_no, 1);
DRO.stage1.lb = zeros(processing_no + edge_no, 1);
DRO.stage1.ub = [max_process_capa_vec; full(capacity_const_vec ...
    - capacity_link_mat * ones(edge_scp_no, 1))];
DRO.stage2 = struct;
DRO.stage2.V = [[sparse(processing_no + supplier_no + customer_no, ...
    processing_no); speye(processing_no); 
    sparse(edge_no, processing_no)], sparse(K2ast, edge_no)];
DRO.stage2.W = blkdiag([sparse(processing_no + supplier_no, ...
    customer_no); speye(customer_no); sparse(processing_no, ...
    customer_no)], -capacity_link_mat);
DRO.stage2.b = [zeros(processing_no, 1); supply_vec; ...
    zeros(customer_no + processing_no, 1); capacity_const_vec];
DRO.stage2.Lin = [p_in_out_mat', s_out_mat', c_in_mat', p_in_mat', ...
    speye(edge_no)];
DRO.stage2.qin = [edges.s2p.cost; edges.p2c.cost];
DRO.stage2.lb = [balance_dual_lb; supply_dual_lb; zeros(customer_no, 1); 
    processing_dual_lb; capacity_dual_lb];
DRO.stage2.ub = [balance_dual_ub; zeros(supplier_no, 1); demand_dual_ub; 
    zeros(processing_no, 1); zeros(edge_no, 1)];

end

