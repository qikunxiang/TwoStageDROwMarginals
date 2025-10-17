classdef (Abstract) TSDROMProb < handle
    % Abstract class for probability measures supported on a one-dimensional interval
    
    properties(Constant)
        % numerical tolerance for deciding if a point is inside an interval
        INSIDE_TOLERANCE = 1e-14;
    end

    properties(SetAccess = protected, GetAccess = public)
        % struct storing information about the support of the probability measure
        Supp = struct;
    end
    
    methods(Access = public)
        function obj = TSDROMProb(intv_lb, intv_ub)
            % Constructor which sets the supporting interval
            % Inputs:
            %   intv_lb: the lower boundary of the interval
            %   intv_ub: the upper boundary of the interval

            assert(intv_ub > intv_lb, 'the interval is empty');

            obj.Supp.LowerBound = intv_lb;
            obj.Supp.UpperBound = intv_ub;
        end
        
        function inside = checkIfInsideInterval(obj, pts)
            % Check if the input points are inside the interval
            % Inputs:
            %   pts: vector containing the input points
            % Output: 
            %   inside: boolean vector indicating whether each point is inside the interval

            inside = pts >= obj.Supp.LowerBound - TSDROMProb.INSIDE_TOLERANCE ...
                & pts <= obj.Supp.UpperBound + TSDROMProb.INSIDE_TOLERANCE;
        end

        function val = computeMean(obj)
            % Compute the mean of the distribution
            % Output:
            %   val: the computed mean value

            val = obj.computeAffIntegral(1, 0, obj.Supp.LowerBound, obj.Supp.UpperBound);
        end

        function samps = randSample(obj, samp_num, rand_stream, varargin)
            % Rejection sampling the probability measure.
            % Inputs:
            %   samp_num: number of samples to genenerate
            %   rand_stream: RandStream object (default is RandStream.getGlobalStream) 
            % Output:
            %   samps: vector containing the generated samples

            if ~exist('rand_stream', 'var') || isempty(rand_stream)
                rand_stream = RandStream.getGlobalStream;
            end

            % first generate samples from Uniform[0, 1]
            u = rand(rand_stream, samp_num, 1);

            % then apply the inverse CDF transform
            samps = obj.computeLContInvCDF(u, varargin{:});
        end

        function samps = randSampleConditional(obj, disc_atoms, disc_probs, disc_samp_indices, rand_stream, varargin)
            % Randomly generate samples from the optimal coupling of the probability measure and a given discrete measure conditional
            % on a set of pre-generated samples from the discrete measure
            % Inputs: 
            %   disc_atoms: vector containing atoms in the discrete measure; assumed to be already sorted into ascending order
            %   disc_probs: vector containing probabilities of the atoms in the discrete measure
            %   disc_samp_indices: vector containing the indices corresponding to the set of pre-generated samples from the discrete
            %   measure
            %   rand_stream: RandStream object (default is RandStream.getGlobalStream) 
            % Output:
            %   samps: vector containing the generated samples
            %   disc_samps: vector containing the corresponding sampled atom indices in the discrete measure

            if ~exist('rand_stream', 'var') || isempty(rand_stream)
                rand_stream = RandStream.getGlobalStream;
            end

            samp_num = length(disc_samp_indices);
            atom_num = length(disc_atoms);
            assert(length(disc_probs) == atom_num, 'discrete measure mis-specified');

            % normalize the discrete measure to prevent numerical error
            disc_probs = disc_probs / sum(disc_probs);

            % compute the cumulative probabilities of atoms from left to right
            cum_probs = cumsum(disc_probs);

            % count the number of samples coupled with each atom
            samp_num_list = accumarray(disc_samp_indices, ones(samp_num, 1), [atom_num, 1]);

            % generate the samples coupled with each atom
            samp_cell = cell(atom_num, 1);

            for atom_id = 1:atom_num
                % the distribution of u after marginalizing out the atom labels in the discrete measure is uniform on [0,1]
                u = cum_probs(atom_id) - disc_probs(atom_id) + disc_probs(atom_id) * rand(rand_stream, samp_num_list(atom_id), 1);
                samp_cell{atom_id} = obj.computeLContInvCDF(u, varargin{:});
            end
            
            samps = vertcat(samp_cell{:});

            % order the samples according the discrete samples
            [~, ordering] = sort(disc_samp_indices, 'ascend');
            samps(ordering) = samps;
        end

        function [OTcost, Kpot_integral, Kpot_dual] = computeOT(obj, disc_atoms, disc_probs)
            % Compute the optimal transport cost to a discrete measure with respect to the cost \R^2 \ni (x, y) \mapsto -x * y \in\R.
            % The Kantorovich potential functions are given by \psi and f, where \psi is a function on the underlying space of this
            % probability measure and f is a function on the underlying space of the discrete measure. In particular, \psi and f 
            % satisfy: \psi(x) + f(z) \ge x * z. Moreover, we force \psi to be 0 when evaluated at the left boundary.
            % Inputs:
            %   disc_atoms: vector containing atoms in the discrete measure; assumed to be already sorted into ascending order
            %   disc_probs: vector containing probabilities of the atoms in the discrete measure
            % Outputs:
            %   OTcost: the computed optimal transport cost
            %   Kpot_integral: the integral of the Kantorovich potential function \psi with respect to this probability measure
            %   Kpot_dual: matrix with four columns where the first column contains the knots in the continuous piece-wise affine
            %   Kantorovich potential function f, the second column contains the values of f evaluated at the knots, and the third and
            %   fourth column contain the left- and right-derivatives of f evaluated at the knots

            atom_num = length(disc_atoms);
            assert(length(disc_probs) == atom_num, 'discrete measure is mis-specified');

            % compute the cumulative probabilities
            cum_probs = cumsum(disc_probs);

            % compute the left-continuous generalized inverse of the CDF; these values serve as subgradients of the function f
            invcdf = obj.computeLContInvCDF(cum_probs(1:end - 1, :));

            % compute the Kantorovich potential function f evaluated at the atoms via the subgradients computed above
            Kpot_dual_vals = cumsum([obj.Supp.LowerBound; invcdf] .* [disc_atoms(1); diff(disc_atoms)]);
            Kpot_dual = [disc_atoms, Kpot_dual_vals, [obj.Supp.LowerBound; invcdf], [invcdf; obj.Supp.UpperBound]];

            % compute the integrals of the Kantorovich potential function \psi over subintervals
            % the right end points of all subintervals but the last one are shifted slightly to the left in order to approximate
            % right-open intervals; this prevents overcounting when there is an atom in the probability measure exactly on a right end
            % point
            Kpot_integral_sub = obj.computeAffIntegral(disc_atoms, -Kpot_dual_vals, ...
                [obj.Supp.LowerBound; invcdf], [invcdf; obj.Supp.UpperBound + 1e-14]);

            % sum up the integral over subintervals
            Kpot_integral = sum(Kpot_integral_sub);

            % the optimal transport cost is given by the negative of the sum of the integral of \psi over this measure and the integral
            % of f over the discrete measure
            OTcost = -Kpot_integral - Kpot_dual_vals' * disc_probs;
        end
    end

    methods(Abstract, Access = public)
        % Compute the left-continuous generalized inverse cumulative density function
        x = computeLContInvCDF(obj, p, varargin);

        % Compute the right limit of the left-continuous generalized inverse cumulative density function
        x = computeLContInvCDFRLim(obj, p, varargin);

        % Compute the integral of an affine function on an interval: \R \ni x \mapsto (a * x + b) * \INDI_{[lb, ub)} \in \R
        vals = computeAffIntegral(obj, a, b, lb, ub);
    end
end

