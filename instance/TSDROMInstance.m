classdef  TSDROMInstance < handle
    %TSDROMINSTANCE Class for instances of the two-stage distributionally robust optimization (DRO) problem with marginal
    %constraints
    
    properties(SetAccess = protected, GetAccess = public)

        % Number of uncertain variables
        unc_num;

        % The object containing the uncertain variables and the marginal constraints
        uncertain;


        % First-stage problem

        % Number of first-stage decision variables 
        st1_deci_num;

        % Objective (cost) vector of the first-stage problem
        st1_obj;

        % Number of inequality constraints in the first-stage problem
        st1_constr_num_in;

        % Inequality constraint matrix of the first-stage problem
        st1_constr_mat_in;

        % Inequality constraint RHS of the first-stage problem
        st1_constr_rhs_in;

        % Number of equality constraints in the first-stage problem
        st1_constr_num_eq;

        % Equality constraint matrix of the first-stage problem
        st1_constr_mat_eq;

        % Equality constraint RHS of the first-stage problem
        st1_constr_rhs_eq;

        % Second-stage problem

        % Number of second-stage decision variables
        st2_deci_num;

        % Objective (cost) vector of the second-stage problem
        st2_obj;

        % Number of inequality constraints in the second-stage problem
        st2_constr_num_in;

        % Inequality constraint matrix of the second-stage problem
        st2_constr_mat_in;

        % Part of the inequality constraint RHS of the second-stage problem that depends on the action (first-stage decision)
        st2_constr_rhs_act_in;

        % Part of the inequality constraint RHS of the second-stage problem that depends on the uncertain variale (random event)
        st2_constr_rhs_unc_in;

        % Part of the inequality constraint RHS of the second-stage problem that is an intercept
        st2_constr_rhs_itc_in;

        % Number of equality constraints in the second-stage problem
        st2_constr_num_eq;

        % Equality constraint matrix of the second-stage problem
        st2_constr_mat_eq;

        % Part of the equality constraint RHS of the second-stage problem that depends on the action (first-stage decision)
        st2_constr_rhs_act_eq;

        % Part of the equality constraint RHS of the second-stage problem that depends on the uncertainty
        st2_constr_rhs_unc_eq;

        % Part of the equality constraint RHS of the second-stage problem that is an intercept
        st2_constr_rhs_itc_eq;

        % Dualized second-stage problem (which is a maximization problem)

        % Number of decision variables in the dualized second-stage problem
        st2d_deci_num;

        % Part of the objective of the dualized second-stage problem that depends on the action (first-stage decision)
        st2d_obj_act;

        % Part of the objective of the dualized second-stage problem that depends on the uncertainty
        st2d_obj_unc;

        % Part of the objective of the dualized second-stage problem that is an intercept
        st2d_obj_itc;

        % Lower bounds for the decision variables in the dualized second-stage problem
        st2d_constr_lb;

        % Upper bounds for the decision variables in the dualized second-stage problem
        st2d_constr_ub;

        % Number of equality constraints in the dualized second-stage problem
        st2d_constr_num_in;

        % Equality constraint matrix of the dualized second-stage problem
        st2d_constr_mat_in;

        % Equality constraint RHS of the dualized second-stage problem
        st2d_constr_rhs_in;

        % Lower bounds for the projections of the decision variables in the dualized second-stage problem by st2d_obj_unc
        st2d_proj_lb;

        % Upper bounds for the projections of the decision variables in the dualized second-stage problem by st2d_obj_unc
        st2d_proj_ub;
    end
    
    methods(Access = public)
        function obj = TSDROMInstance( ...
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
                st2_constr_rhs_itc_eq)
            % Constructor

            obj.unc_num = uncertain.num;
            obj.uncertain = uncertain;
            obj.uncertain.setProblem(obj);

            % check for dimensionality of inputs
            obj.st1_deci_num = st1_deci_num;
            assert(size(st1_obj, 1) == obj.st1_deci_num && size(st1_obj, 2) == 1, ...
                'st1_obj invalid');
            obj.st1_obj = st1_obj;
            
            obj.st1_constr_num_in = st1_constr_num_in;
            assert(size(st1_constr_mat_in, 1) == obj.st1_constr_num_in && size(st1_constr_mat_in, 2) == obj.st1_deci_num, ...
                'st1_constr_mat_in invalid');
            obj.st1_constr_mat_in = st1_constr_mat_in;
            assert(size(st1_constr_rhs_in, 1) == obj.st1_constr_num_in && size(st1_constr_rhs_in, 2) == 1, ...
                'st1_constr_rhs_in invalid');
            obj.st1_constr_rhs_in = st1_constr_rhs_in;

            obj.st1_constr_num_eq = st1_constr_num_eq;
            assert(size(st1_constr_mat_eq, 1) == obj.st1_constr_num_eq && size(st1_constr_mat_eq, 2) == obj.st1_deci_num, ...
                'st1_constr_mat_eq invalid');
            obj.st1_constr_mat_eq = st1_constr_mat_eq;
            assert(size(st1_constr_rhs_eq, 1) == obj.st1_constr_num_eq && size(st1_constr_rhs_eq, 2) == 1, ...
                'st1_constr_rhs_eq invalid');
            obj.st1_constr_rhs_eq = st1_constr_rhs_eq;

            obj.st2_deci_num = st2_deci_num;
            assert(size(st2_obj, 1) == obj.st2_deci_num && size(st2_obj, 2) == 1, ...
                'st2_obj invalid');
            obj.st2_obj = st2_obj;

            obj.st2_constr_num_in = st2_constr_num_in;
            assert(size(st2_constr_mat_in, 1) == obj.st2_constr_num_in && size(st2_constr_mat_in, 2) == obj.st2_deci_num, ...
                'st2_constr_mat_in invalid');
            obj.st2_constr_mat_in = st2_constr_mat_in;
            assert(size(st2_constr_rhs_act_in, 1) == obj.st2_constr_num_in && size(st2_constr_rhs_act_in, 2) == obj.st1_deci_num, ...
                'st2_constr_rhs_act_in invalid');
            obj.st2_constr_rhs_act_in = st2_constr_rhs_act_in;
            assert(size(st2_constr_rhs_unc_in, 1) == obj.st2_constr_num_in && size(st2_constr_rhs_unc_in, 2) == obj.unc_num, ...
                'st2_constr_rhs_unc_in invalid');
            obj.st2_constr_rhs_unc_in = st2_constr_rhs_unc_in;
            assert(size(st2_constr_rhs_itc_in, 1) == obj.st2_constr_num_in && size(st2_constr_rhs_itc_in, 2) == 1, ...
                'st2_constr_rhs_itc_in invalid');
            obj.st2_constr_rhs_itc_in = st2_constr_rhs_itc_in;

            obj.st2_constr_num_eq = st2_constr_num_eq;
            assert(size(st2_constr_mat_eq, 1) == obj.st2_constr_num_eq && size(st2_constr_mat_eq, 2) == obj.st2_deci_num, ...
                'st2_constr_mat_eq invalid');
            obj.st2_constr_mat_eq = st2_constr_mat_eq;
            assert(size(st2_constr_rhs_act_eq, 1) == obj.st2_constr_num_eq && size(st2_constr_rhs_act_eq, 2) == obj.st1_deci_num, ...
                'st2_constr_rhs_act_eq invalid');
            obj.st2_constr_rhs_act_eq = st2_constr_rhs_act_eq;
            assert(size(st2_constr_rhs_unc_eq, 1) == obj.st2_constr_num_eq && size(st2_constr_rhs_unc_eq, 2) == obj.unc_num, ...
                'st2_constr_rhs_unc_eq invalid');
            obj.st2_constr_rhs_unc_eq = st2_constr_rhs_unc_eq;
            assert(size(st2_constr_rhs_itc_eq, 1) == obj.st2_constr_num_eq && size(st2_constr_rhs_itc_eq, 2) == 1, ...
                'st2_constr_rhs_itc_eq invalid');
            obj.st2_constr_rhs_itc_eq = st2_constr_rhs_itc_eq;

            
            % Initialize the dualized second-stage problem
            obj.st2d_deci_num = obj.st2_constr_num_in + obj.st2_constr_num_eq;
            obj.st2d_obj_act = [obj.st2_constr_rhs_act_in; obj.st2_constr_rhs_act_eq];
            obj.st2d_obj_unc = [obj.st2_constr_rhs_unc_in; obj.st2_constr_rhs_unc_eq];
            obj.st2d_obj_itc = [obj.st2_constr_rhs_itc_in; obj.st2_constr_rhs_itc_eq];
            obj.st2d_constr_lb = -inf(obj.st2d_deci_num, 1);
            obj.st2d_constr_ub = [zeros(obj.st2_constr_num_in, 1); inf(obj.st2_constr_num_eq, 1)];
            obj.st2d_constr_num_in = obj.st2_deci_num;
            obj.st2d_constr_mat_in = [obj.st2_constr_mat_in' obj.st2_constr_mat_eq'];
            obj.st2d_constr_rhs_in = obj.st2_obj;

            % by default there are no known lower and upper bounds for the projection of the decision variable of the dualized
            % second-stage problem
            obj.st2d_proj_lb = -inf(obj.unc_num, 1);
            obj.st2d_proj_ub = inf(obj.unc_num, 1);
        end
    end
end

