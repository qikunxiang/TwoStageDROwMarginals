classdef TSDROMSupplyNetUncertain < TSDROMUncertain
    % Class representing the uncertain parameters in the two-stage distributionally robust supply chain network design problem with 
    % marginal constraints
    
    properties(SetAccess = protected, GetAccess = public)
        demand_marginals;

        edge_failure_marginals;

        demand_max = [];

        customer_num;

        edge_susceptible_num;

        edge_failure_max;
    end
    
    methods(Access = public)
        function obj = TSDROMSupplyNetUncertain(demand_marginals, edge_failure_marginals)
            % Constructor
            % Inputs:
            %   demand_marginals: cell array containing the marginals representing the demands
            %   edge_failure_marginals: cell array containing the marginals representing the failure of edges

            customer_num = length(demand_marginals);
            edge_susceptible_num = length(edge_failure_marginals);

            demand_max = zeros(customer_num, 1);

            for cust_id = 1:customer_num
                demand_max(cust_id) = demand_marginals{cust_id}.Supp.UpperBound;
            end

            edge_failure_max = zeros(edge_susceptible_num, 1);

            for edge_id = 1:edge_susceptible_num
                marg = edge_failure_marginals{edge_id};
                assert(marg.Supp.LowerBound == 0 && marg.Supp.UpperBound <= 1, 'each edge failure should be between 0 and 1');
                edge_failure_max(edge_id) = marg.Supp.UpperBound;
            end

            obj@TSDROMUncertain([demand_marginals; edge_failure_marginals]);

            obj.demand_marginals = demand_marginals;
            obj.edge_failure_marginals = edge_failure_marginals;
            obj.customer_num = customer_num;
            obj.edge_susceptible_num = edge_susceptible_num;
            obj.demand_max = demand_max;
            obj.edge_failure_max = edge_failure_max;
        end
    end
end

