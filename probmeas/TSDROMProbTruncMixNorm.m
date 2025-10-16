classdef TSDROMProbTruncMixNorm < TSDROMProb
    % Class for probability measures supported on a one-dimensional interval which is a mixture of normal distributions truncated to
    % the interval
    
    properties(SetAccess = protected, GetAccess = public)
        % struct containing the parameters
        Params = struct;

        % struct containing options
        Options = struct;
    end
    
    methods(Access = public)
        function obj = TSDROMProbTruncMixNorm(intv_lb, intv_ub, weight_list, mean_list, std_list, opt)
            % Constructor
            % Inputs:
            %   intv_lb: left boundary of the support interval
            %   intv_ub: right boundary of the support interval
            %   weight_list: weights of the mixture components
            %   mean_list: vector containing the mean parameters of the mixture components
            %   std_list: vector containing the standard deviation parameters of the mixture components
            %   opt: struct containing the options

            obj@TSDROMProb(intv_lb, intv_ub);

            obj.Params.comp_num = size(weight_list, 1);

            assert(size(weight_list, 2) == 1, 'weights mis-specified');

            obj.Params.weight_list = weight_list;

            assert(size(mean_list, 1) == obj.Params.comp_num && size(mean_list, 2) == 1, 'mean_list mis-specified');
            assert(size(std_list, 1) == obj.Params.comp_num && size(std_list, 2) == 1, 'std_list mis-specified');

            assert(all(std_list >= 1e-6 | (mean_list >= intv_lb & mean_list <= intv_ub)), ...
                'there are degenerate components outside the support inteval');

            if all(std_list < 1e-6)
                % in the case where all components are degenerate, shrink the support interval
                obj.Supp.LowerBound = min(mean_list);
                obj.Supp.UpperBound = max(mean_list);
            end

            obj.Params.mean_list = mean_list;
            obj.Params.std_list = std_list;

            obj.Params.untruncated_cdf_lb = obj.computeUntruncatedCDF(intv_lb);
            obj.Params.untruncated_cdf_ub = obj.computeUntruncatedCDF(intv_ub);
            obj.Params.normalize_const = obj.Params.untruncated_cdf_ub - obj.Params.untruncated_cdf_lb;

            if ~exist('opt', 'var') || isempty(opt)
                opt = struct;
            end

            obj.Options = opt;

            if ~isfield(obj.Options, 'bisect_bin_num') || isempty(obj.Options.bisect_bin_num)
                % number of bins used to initialize the bisection method
                obj.Options.bisect_bin_num = 2^10;
            end

            if ~isfield(obj.Options, 'bisect_iter_num') || isempty(obj.Options.bisect_iter_num)
                % number of iterations in the bisection method
                % for example, when obj.Options.bisect_bin_num = 2^10 and obj.Options.bisect_iter_num = 16, the maximum error in the
                % output of the bisection method is equal to 2^-17 = 5.4506e-9
                obj.Options.bisect_iter_num = 16;
            end

            if ~isfield(obj.Options, 'boundary_prob_tolerance') || isempty(obj.Options.boundary_prob_tolerance)
                % when evaluating the inverse CDF with inputs close to 0/1, the returned values might not be the lower/upper bound due 
                % to the probabilities being "concentrated" away from the boundary; to resolve this, we impose that the inverse CDF of
                % probabilities whose distances to 0/1 fall below this threshold to be equal to the corresponding boundaries
                obj.Options.boundary_prob_tolerance = 1e-15;
            end

            obj.initializeBisection();
        end

        function p = computeUntruncatedCDF(obj, x)
            % Compute the cumulative density function without truncation
            % Input:
            %   x: vector containing the inputs
            % Output:
            %   p: vector containing the computed probabilities

            p = sum((normcdf((x - obj.Params.mean_list') ./ obj.Params.std_list')) .* obj.Params.weight_list', 2);
        end

        function p = computeCDF(obj, x)
            % Compute the cumulative density function
            % Input:
            %   x: vector containing the inputs
            % Output:
            %   p: vector containing the computed probabilities

            p = min(max((obj.computeUntruncatedCDF(x) - obj.Params.untruncated_cdf_lb) / obj.Params.normalize_const, 0), 1);
        end
        
        function x = computeLContInvCDF(obj, p)
            % Compute the left-continuous generalized inverse cumulative density function
            % Input:
            %   p: vector containing the probabilities with respect to which the inverse CDF is to be evaluated
            % Output:
            %   x: vector containing the computed inverse CDF

            p = min(max(p, 0), 1);

            % compute the corresponding probabilities in the untruncated distribution
            p_ut = p * obj.Params.normalize_const + obj.Params.untruncated_cdf_lb;

            if obj.Params.comp_num == 1
                % if there is only a single component, then norminv can be used
                x = norminv(p_ut * obj.Params.normalize_const + obj.Params.cdf_lb) * obj.Params.std_list + obj.Params.mean_list;
            else
                % else, use the bisection method to compute the inverse CDF

                % initial intervals
                init_bin_indices = discretize(p_ut, obj.Params.bisect_discretized_vals);
                x_ranges = [obj.Params.bisect_bin_edges(init_bin_indices), obj.Params.bisect_bin_edges(init_bin_indices + 1)];

                % perform bisection for a fixed number of times
                for bisect_iter = 1:obj.Options.bisect_iter_num
                    x_mid = mean(x_ranges, 2);
                    p_mid = obj.computeUntruncatedCDF(x_mid);
                    left_list = p_mid <= p_ut;
                    right_list = ~left_list;
                    x_ranges(left_list, 1) = x_mid(left_list);
                    x_ranges(right_list, 2) = x_mid(right_list);
                end

                x = mean(x_ranges, 2);
            end

            x(p < obj.Options.boundary_prob_tolerance) = obj.Supp.LowerBound;
            x(p > 1 - obj.Options.boundary_prob_tolerance) = obj.Supp.UpperBound;
        end

        function x = computeLContInvCDFRLim(obj, p)
            % Compute the right limit of the left-continuous generalized inverse cumulative density function
            % Input:
            %   p: vector containing the probabilities with respect to which the inverse CDF is to be evaluated
            % Output:
            %   x: vector containing the computed inverse CDF

            x = obj.computeLContInvCDF(min(p + 1e-6, 1));
        end

        function vals = computeAffIntegral(obj, a, b, lb, ub)
            % Compute the integral of an affine function on an interval: \R \ni x \mapsto (a * x + b) * \INDI_{[lb, ub]} \in \R
            % Inputs:
            %   a, b, lb, ub: vectors of the same length representing the parameters in the integrand

            term_num = size(a, 1);
            assert(size(a, 2) == 1 && size(b, 1) == term_num && size(b, 2) == 1 && size(lb, 1) == term_num && size(lb, 2) == 1 ...
                && size(ub, 1) == term_num && size(ub, 2) == 1, 'inputs mis-specified');

            vals = zeros(term_num, 1);

            lb = max(lb, obj.Supp.LowerBound);
            ub = max(min(ub, obj.Supp.UpperBound), lb);

            for comp_id = 1:obj.Params.comp_num
                comp_mean = obj.Params.mean_list(comp_id);
                comp_std = obj.Params.std_list(comp_id);

                % standardize the lower and upper integration limits
                lb_standard = (lb - comp_mean) / comp_std;
                ub_standard = (ub - comp_mean) / comp_std;
                
                diff_normcdf = normcdf(ub_standard) - normcdf(lb_standard);
                diff_exp = exp(-0.5 * (lb_standard .^ 2)) - exp(-0.5 * (ub_standard .^ 2));

                % for degenerate components, diff_exp will evaluate to inf, -inf, or nan; we set them to 0
                diff_exp(isinf(diff_exp) | isnan(diff_exp)) = 0;

                norm_partialexp = (a * comp_mean + b) .* diff_normcdf + (comp_std / sqrt(2 * pi)) * a .* diff_exp;

                vals = vals + norm_partialexp * obj.Params.weight_list(comp_id);
            end

            vals = vals / obj.Params.normalize_const;
        end
    end

    methods(Access = protected)
        function initializeBisection(obj)
            % Initialize the precomputed quantities used in the bisection method

            obj.Params.bisect_bin_edges = (obj.Supp.UpperBound - obj.Supp.LowerBound) / obj.Options.bisect_bin_num ...
                * (0:obj.Options.bisect_bin_num)' + obj.Supp.LowerBound;
            obj.Params.bisect_discretized_vals = obj.computeUntruncatedCDF(obj.Params.bisect_bin_edges);
            obj.Params.bisect_discretized_vals(1) = obj.Params.untruncated_cdf_lb - eps;
            obj.Params.bisect_discretized_vals(end) = obj.Params.untruncated_cdf_ub + eps;
        end
    end
end

