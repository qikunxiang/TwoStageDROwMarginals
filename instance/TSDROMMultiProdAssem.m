classdef TSDROMMultiProdAssem < TSDROMInstance
    % Class for the two-stage distributionally robust multi-product assembly problem with marginal constraints
    
    properties(SetAccess = protected, GetAccess = public)
        % Number of parts
        part_num;

        % Number of products
        product_num;

        % Vector containing the costs of parts
        part_costs;

        % Matrix describing the parts needed for producing each unit of product; each row represents a product and each column 
        % represents a part
        assembly_parts;

        % Vector containing the revenues from selling one unit of each type of product
        product_revenues;

        % Vector containing the salvage prices of the remaining parts
        part_salvage_prices;
    end
    
    methods(Access = public)
        function obj = TSDROMMultiProdAssem(uncertain, part_costs, assembly_parts, product_revenues, part_salvage_prices)
            % Constructor
            % Inputs: 
            %   part_costs: vector containing the costs of parts
            %   assembly_parts: matrix describing the parts needed for producing each unit of product; each row represents a product 
            %   and each column represents a part
            %   product_revenues: vector containing the revenues from selling one unit of each type of product
            %   part_salvage_prices: vector containing the salvage prices of the remaining parts

            product_num = uncertain.num;
            part_num = size(part_costs, 1);
            assert(size(part_costs, 2) == 1, 'part_costs mis-specified');
            assert(size(assembly_parts, 1) == product_num && size(assembly_parts, 2) == part_num, 'assembly_parts mis-specified');
            assert(size(product_revenues, 1) == product_num && size(product_revenues, 2) == 1, 'product_revenues mis-specified');
            assert(size(part_salvage_prices, 1) == part_num && size(part_salvage_prices, 2) == 1, 'part_salvage_prices mis-specified');

            assembly_parts = sparse(assembly_parts);

            st1_deci_num = part_num;
            st1_obj = part_costs;
            st1_constr_num_in = 0;
            st1_constr_mat_in = speye(0, part_num);
            st1_constr_rhs_in = zeros(0, 1);
            st1_constr_num_eq = 0;
            st1_constr_mat_eq = sparse(0, part_num);
            st1_constr_rhs_eq = zeros(0, 1);
            st2_deci_num = part_num + product_num;
            st2_obj = [-product_revenues; -part_salvage_prices];
            st2_constr_num_in = product_num;
            st2_constr_mat_in = [speye(product_num), sparse(product_num, part_num)];
            st2_constr_rhs_act_in = sparse(st2_constr_num_in, st1_deci_num);
            st2_constr_rhs_unc_in = speye(product_num);
            st2_constr_rhs_itc_in = zeros(st2_constr_num_in, 1);
            st2_constr_num_eq = part_num;
            st2_constr_mat_eq = [assembly_parts', speye(part_num)];
            st2_constr_rhs_act_eq = speye(part_num);
            st2_constr_rhs_unc_eq = sparse(part_num, product_num);
            st2_constr_rhs_itc_eq = zeros(part_num, 1);

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

            obj.part_num = part_num;
            obj.product_num = product_num;
            obj.part_costs = part_costs;
            obj.assembly_parts = assembly_parts;
            obj.product_revenues = product_revenues;
            obj.part_salvage_prices = part_salvage_prices;
            obj.st2d_proj_ub = zeros(obj.unc_num, 1);
        end
    end
end

