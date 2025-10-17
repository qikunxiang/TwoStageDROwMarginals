classdef TSDROMProbBernoulli < TSDROMProb
    % Class for probability measures corresponding to Bernoulli distributions
    
    properties(GetAccess = public, SetAccess = protected)
        probability = struct;
    end
    
    methods(Access = public)
        function obj = TSDROMProbBernoulli(prob)
            % Constructor
            %   Inputs:
            %       prob: probability of 1
            
            assert(prob >= 0 && prob <= 1, 'probability must be between 0 and 1');

            obj@TSDROMProb(0, 1);
            obj.probability = prob;
        end

        function x = computeLContInvCDF(obj, p)
            % Compute the left-continuous generalized inverse cumulative density function
            % Input:
            %   p: vector containing the probabilities with respect to which the inverse CDF is to be evaluated
            % Output:
            %   x: vector containing the computed inverse CDF

            x = double(p > (1 - obj.probability));
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

            vals = sum((a .* [0, 1] + b) .* (lb <= [0, 1] & ub > [0, 1]) .* [1 - obj.probability, obj.probability], 2);
        end
    end
end

