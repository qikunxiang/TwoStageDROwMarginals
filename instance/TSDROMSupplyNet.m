classdef TSDROMSupplyNet < TSDROMInstance
    % Class for the two-stage distributionally robust supply chain network design problem with marginal constraints
    
    properties(SetAccess = protected, GetAccess = public)
        % Number of suppliers
        supplier_num;

        % Number of processing facilities
        process_num;

        % Number of customers
        customer_num;

        % Number of edges from the suppliers to the processing facilities
        edge_S2P_num;

        % Number of edges from the processing facilities to the customers
        edge_P2C_num;

        % Number of edges in total
        edge_num;

        % Nubmer of edges susceptible to failure
        edge_susceptible_num;

        % Struct representing the network
        network;

        % Supplies of the suppliers
        supplies;

        % Maximum processing capabilities of the processing facilities
        process_capabilities_max;

        % Per unit cost of investment into the processing capabilities of the processing facilities
        process_investment_costs;

        % Maximum demand of the customers
        demand_max;

        % Minimum capacities of the failed edges
        edge_failure_max;
    end
    
    methods(Access = public)
        function obj = TSDROMSupplyNet(uncertain, ...
                network, ...
                supplies, ...
                process_capabilities_max, ...
                process_investment_costs)
            % Constructor
            % Inputs: 
            %   uncertain: instance of TsDROMSupplyNetUncertain
            %   network: struct containing the following fields
            %       S2P: struct containing the fields sources, destinations, costs, capacities, and susceptibilities that are vectors 
            %       of equal lengths
            %       P2C: struct containing the fields sources, destinations, costs, capacities, and susceptibilities that are vectors 
            %       of equal lengths
            %   supplies: vector containing the supply of each supplier (sum of the supply must exceed the sum of the maximum demand)
            %   process_capabilities_max: vector containing the maximum processing capabilities of the processing facilities
            %   process_investment_costs: vector containing the per unit cost of investment into the processing capability of each 
            %   processing facility

            assert(isa(uncertain, 'TSDROMSupplyNetUncertain'), 'uncertain needs to be an instance of TSDROMSupplyNetUncertain');

            supplier_num = size(supplies, 1);
            assert(size(supplies, 2) == 1, 'supplies mis-specified');
            process_num = size(process_capabilities_max, 1);
            assert(size(process_capabilities_max, 2) == 1, 'process_capabilities_max mis-specified');
            assert(size(process_investment_costs, 1) == process_num && size(process_investment_costs, 2) == 1, ...
                'process_investment_costs mis-specified');
            customer_num = uncertain.customer_num;
            demand_max = uncertain.demand_max;
            edge_failure_max = uncertain.edge_failure_max;

            assert(sum(supplies) >= sum(demand_max), 'total supply should not be less than the maximum possible demand');
            assert(sum(process_capabilities_max) >= sum(demand_max), ...
                'maximum processing capability should not be less than the maximum possible demand');

            edge_S2P_num = size(network.S2P.sources, 1);
            assert(size(network.S2P.sources, 2) == 1 ...
                && size(network.S2P.destinations, 1) == edge_S2P_num && size(network.S2P.destinations, 2) == 1 ...
                && min(network.S2P.sources) >= 1 && max(network.S2P.sources) <= supplier_num ...
                && min(network.S2P.destinations) >= 1 && max(network.S2P.destinations) <= process_num ...
                && size(network.S2P.costs, 1) == edge_S2P_num && size(network.S2P.costs, 2) == 1 ...
                && size(network.S2P.capacities, 1) == edge_S2P_num && size(network.S2P.capacities, 2) == 1 ...
                && size(network.S2P.susceptibilities, 1) == edge_S2P_num && size(network.S2P.susceptibilities, 2) == 1, ...
                'S2P edges in the network mis-specified');
            edge_P2C_num = size(network.P2C.costs, 1);
            assert(size(network.P2C.sources, 2) == 1 ...
                && size(network.P2C.destinations, 1) == edge_P2C_num && size(network.P2C.destinations, 2) == 1 ...
                && min(network.P2C.sources) >= 1 && max(network.P2C.sources) <= process_num ...
                && min(network.P2C.destinations) >= 1 && max(network.P2C.destinations) <= customer_num ...
                && size(network.P2C.costs, 1) == edge_P2C_num && size(network.P2C.costs, 2) == 1 ...
                && size(network.P2C.capacities, 1) == edge_P2C_num && size(network.P2C.capacities, 2) == 1 ...
                && size(network.P2C.susceptibilities, 1) == edge_P2C_num && size(network.P2C.susceptibilities, 2) == 1, ...
                'P2C edges in the network mis-specified');
            edge_num = edge_S2P_num + edge_P2C_num;
            edge_susceptible_num = sum(network.S2P.susceptibilities) + sum(network.P2C.susceptibilities);
            assert(uncertain.edge_susceptible_num == edge_susceptible_num, 'uncertain and network do not match');

            % sparse binary matrix where each row represents a supplier and each entry represents an edge exiting the supplier
            S2P_mat = sparse(network.S2P.sources, (1:edge_S2P_num)', 1, supplier_num, edge_S2P_num);

            % sparse binary matrix where each row represents a processing facility and each entry represents an edge entering the
            % processing facility
            P2S_mat = sparse(network.S2P.destinations, (1:edge_S2P_num)', 1, process_num, edge_S2P_num);

            % sparse binary matrix where each row represents a processing facility and each entry represents an edge exiting the
            % processing facility
            P2C_mat = sparse(network.P2C.sources, (1:edge_P2C_num)', 1, process_num, edge_P2C_num);

            % sparse binary matrix where each row represents a customer and each entry represents an edge entering the customer
            C2P_mat = sparse(network.P2C.destinations, (1:edge_P2C_num)', 1, customer_num, edge_P2C_num);

            % sparse matrix representing the balance equality constraints
            balance_constr_mat = [P2S_mat, -P2C_mat];

            % sparse matrix representing the supply inequality constraints
            supply_constr_mat = [S2P_mat, sparse(supplier_num, edge_P2C_num)];

            % sparse matrix representing the demand inequality constraints
            demand_constr_mat = [sparse(customer_num, edge_S2P_num), C2P_mat];

            % sparse matrix representing the processing capability inequality constraints
            process_constr_mat = [P2S_mat, sparse(process_num, edge_P2C_num)];

            % sparse matrix representing the susceptible edges and the corresponding capacities
            edge_capacities = [network.S2P.capacities; network.P2C.capacities];
            edge_susceptible_indices = find([network.S2P.susceptibilities; network.P2C.susceptibilities]);
            susceptible_mat = sparse(edge_susceptible_indices, (1:edge_susceptible_num)', ...
                edge_capacities(edge_susceptible_indices), edge_num, edge_susceptible_num);

            st1_deci_num = process_num + edge_num;
            st1_obj = [process_investment_costs; zeros(edge_num, 1)];

            % stage 1 constraints consist of 6 parts:
            %   i. the processing capabilities are bounded above by the maximum processing capabilities
            %  ii. the auxiliary network flow are bounded above by the edge capacities with maximum failure
            % iii. the auxiliary network flow needs to satisfy the supply inequality constraints
            %  iv. the auxiliary network flow needs to satisfy the demand inequality constraints
            %   v. the processing capabilities and the auxiliary network flow need to satisfy the processing capability inequality
            %   constraints
            %  vi. the auxiliary network flow needs to satisfy the balance equality constraints

            st1_constr_num_in_process_ub = process_num;
            st1_constr_mat_in_process_ub = [speye(process_num), sparse(process_num, edge_num)];
            st1_constr_rhs_in_process_ub = process_capabilities_max;

            st1_constr_num_in_flow_capacity = edge_num;
            st1_constr_mat_in_flow_capacity = [sparse(edge_num, process_num), speye(edge_num)];
            st1_constr_rhs_in_flow_capacity = [network.S2P.capacities; network.P2C.capacities] - susceptible_mat * edge_failure_max;

            st1_constr_num_in_flow_supply = supplier_num;
            st1_constr_mat_in_flow_supply = [sparse(supplier_num, process_num), supply_constr_mat];
            st1_constr_rhs_in_flow_supply = supplies;

            st1_constr_num_in_flow_demand = customer_num;
            st1_constr_mat_in_flow_demand = [sparse(customer_num, process_num), -demand_constr_mat];
            st1_constr_rhs_in_flow_demand = -demand_max;

            st1_constr_num_in_flow_process = process_num;
            st1_constr_mat_in_flow_process = [-speye(process_num), process_constr_mat];
            st1_constr_rhs_in_flow_process = zeros(process_num, 1);

            st1_constr_num_in = st1_constr_num_in_process_ub + st1_constr_num_in_flow_capacity + st1_constr_num_in_flow_supply ...
                + st1_constr_num_in_flow_demand + st1_constr_num_in_flow_process;
            st1_constr_mat_in = [st1_constr_mat_in_process_ub; st1_constr_mat_in_flow_capacity; st1_constr_mat_in_flow_supply; ...
                st1_constr_mat_in_flow_demand; st1_constr_mat_in_flow_process];
            st1_constr_rhs_in = [st1_constr_rhs_in_process_ub; st1_constr_rhs_in_flow_capacity; st1_constr_rhs_in_flow_supply; ...
                st1_constr_rhs_in_flow_demand; st1_constr_rhs_in_flow_process];

            st1_constr_num_eq = process_num;
            st1_constr_mat_eq = [sparse(process_num, process_num), balance_constr_mat];
            st1_constr_rhs_eq = zeros(process_num, 1);


            st2_deci_num = edge_num;
            st2_obj = [network.S2P.costs; network.P2C.costs];

            % stage 2 constraints consist of 5 parts:
            %   i. the network flow is bounded above by the edge capacities
            %  ii. the supply inequality constraints
            % iii. the demand inequality constraints
            %  iv. the processing capability inequality constraints
            %   v. the balance equality constraints

            st2_constr_num_in_capacity = edge_num;
            st2_constr_mat_in_capacity = speye(edge_num);
            st2_constr_rhs_act_in_capacity = sparse(edge_num, st1_deci_num);
            st2_constr_rhs_unc_in_capacity = [sparse(edge_num, customer_num), -susceptible_mat];
            st2_constr_rhs_itc_in_capacity = [network.S2P.capacities; network.P2C.capacities];
            
            st2_constr_num_in_supply = supplier_num;
            st2_constr_mat_in_supply = supply_constr_mat;
            st2_constr_rhs_act_in_supply = sparse(supplier_num, st1_deci_num);
            st2_constr_rhs_unc_in_supply = sparse(supplier_num, customer_num + edge_susceptible_num);
            st2_constr_rhs_itc_in_supply = supplies;

            st2_constr_num_in_demand = customer_num;
            st2_constr_mat_in_demand = -demand_constr_mat;
            st2_constr_rhs_act_in_demand = sparse(customer_num, st1_deci_num);
            st2_constr_rhs_unc_in_demand = [-speye(customer_num), sparse(customer_num, edge_susceptible_num)];
            st2_constr_rhs_itc_in_demand = zeros(customer_num, 1);

            st2_constr_num_in_process = process_num;
            st2_constr_mat_in_process = process_constr_mat;
            st2_constr_rhs_act_in_process = [speye(process_num), sparse(process_num, edge_num)];
            st2_constr_rhs_unc_in_process = sparse(process_num, customer_num + edge_susceptible_num);
            st2_constr_rhs_itc_in_process = zeros(process_num, 1);

            st2_constr_num_in = st2_constr_num_in_capacity + st2_constr_num_in_supply + st2_constr_num_in_demand ...
                + st2_constr_num_in_process;
            st2_constr_mat_in = [st2_constr_mat_in_capacity; st2_constr_mat_in_supply; st2_constr_mat_in_demand; ...
                st2_constr_mat_in_process];
            st2_constr_rhs_act_in = [st2_constr_rhs_act_in_capacity; st2_constr_rhs_act_in_supply; st2_constr_rhs_act_in_demand; ...
                st2_constr_rhs_act_in_process];
            st2_constr_rhs_unc_in = [st2_constr_rhs_unc_in_capacity; st2_constr_rhs_unc_in_supply; st2_constr_rhs_unc_in_demand; ...
                st2_constr_rhs_unc_in_process];
            st2_constr_rhs_itc_in = [st2_constr_rhs_itc_in_capacity; st2_constr_rhs_itc_in_supply; st2_constr_rhs_itc_in_demand; ...
                st2_constr_rhs_itc_in_process];

            st2_constr_num_eq = process_num;
            st2_constr_mat_eq = balance_constr_mat;
            st2_constr_rhs_act_eq = sparse(process_num, st1_deci_num);
            st2_constr_rhs_unc_eq = sparse(process_num, customer_num + edge_susceptible_num);
            st2_constr_rhs_itc_eq = zeros(process_num, 1);

            obj@TSDROMInstance( ...
                uncertain, ...
                st1_deci_num, ...
                st1_obj, ...
                st1_constr_num_in, ...
                st1_constr_mat_in, ...
                st1_constr_rhs_in, ...
                st1_constr_num_eq, ...
                st1_constr_mat_eq, ...
                st1_constr_rhs_eq, ...
                st2_deci_num, ...
                st2_obj, ...
                st2_constr_num_in, ...
                st2_constr_mat_in, ...
                st2_constr_rhs_act_in, ...
                st2_constr_rhs_unc_in, ...
                st2_constr_rhs_itc_in, ...
                st2_constr_num_eq, ...
                st2_constr_mat_eq, ...
                st2_constr_rhs_act_eq, ...
                st2_constr_rhs_unc_eq, ...
                st2_constr_rhs_itc_eq);

            obj.supplier_num = supplier_num;
            obj.process_num = process_num;
            obj.customer_num = customer_num;
            obj.edge_S2P_num = edge_S2P_num;
            obj.edge_P2C_num = edge_P2C_num;
            obj.edge_num = edge_num;
            obj.edge_susceptible_num = edge_susceptible_num;
            obj.network = network;
            obj.supplies = supplies;
            obj.process_capabilities_max = process_capabilities_max;
            obj.process_investment_costs = process_investment_costs;
            obj.demand_max = demand_max;
            obj.edge_failure_max = edge_failure_max;

            obj.st2d_proj_lb = zeros(customer_num + edge_susceptible_num, 1);
        end

        function feas = checkFeasibility(obj)
            % Check whether the problem instance is feasible
            % Output:
            %   feas: boolean indicating feasibility of the problem instance

            model = struct;
            model.modelsense = 'min';
            model.obj = zeros(obj.st1_deci_num, 1);
            model.objcon = 0;
            model.A = [obj.st1_constr_mat_in; obj.st1_constr_mat_eq];
            model.rhs = [obj.st1_constr_rhs_in; obj.st1_constr_rhs_eq];
            model.sense = [repmat('<', obj.st1_constr_num_in, 1); repmat('=', obj.st1_constr_num_eq, 1)];

            params = struct;
            params.OutputFlag = 0;

            result = gurobi(model, params);
            feas = strcmp(result.status, 'OPTIMAL');
        end
    end
end

