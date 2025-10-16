classdef TSDROMWorstCaseSampler < handle
    % Class for generating samples from the computed worst-case probability measure associated with the computed first-stage decision
    
    properties(Constant)
        % Constant added to every element before determining the approximately unique elements
        UNIQUENESS_SHIFT = exp(1);

        % The number of digits to round each element to determine the approximately unique elements
        UNIQUENESS_ROUNDING = 5;
    end

    properties(SetAccess = protected, GetAccess = public)
        % The problem instance of class TSDROMInstance
        problem;

        % First-stage decision
        st1_deci;

        % Probabilities over the dual points
        dual_points_probs;

        % Cell array where each cell is a vector containing the indices mapping the original dual points to the projections
        proj_mappings;

        % Cell array where each cell is a vector containing the projected probabilities
        proj_probs;

        % Cell array where each cell is a vector containing the projected probabilities
        proj_atoms;

        % Random stream; default is RandStream.getGlobalStream
        rand_stream;
    end
    
    methods
        function obj = TSDROMWorstCaseSampler(problem, primal_sol, dual_sol)
            % Constructor
            % Inputs:
            %   problem: the problem instance of class TSDROMInstance
            %   primal_sol: primal solution
            %   dual_sol: dual solution

            obj.problem = problem;
            obj.rand_stream = RandStream.getGlobalStream;

            obj.st1_deci = primal_sol.st1_deci;
            dual_points = dual_sol.dual_points;
            probs = dual_sol.dual_points_probs;

            assert(length(obj.st1_deci) == obj.problem.st1_deci_num, 'dimensionality of first-stage decision is incoorect');
            assert(size(dual_points, 1) == length(probs) && size(dual_points, 2) == obj.problem.st2d_deci_num, ...
                'dimensionality of dual points is incorrect');

            obj.dual_points_probs = probs;

            obj.proj_atoms = cell(obj.problem.unc_num, 1);
            obj.proj_mappings = cell(obj.problem.unc_num, 1);

            for marg_id = 1:obj.problem.unc_num
                proj_i = full(dual_points * obj.problem.st2d_obj_unc(:, marg_id));

                [~, indices, mapping] = unique(round(proj_i + obj.UNIQUENESS_SHIFT, obj.UNIQUENESS_ROUNDING), 'sorted');
                obj.proj_atoms{marg_id} = proj_i(indices);
                obj.proj_mappings{marg_id} = mapping;

                obj.proj_probs{marg_id} = accumarray(obj.proj_mappings{marg_id}, probs, size(obj.proj_atoms{marg_id}));
                obj.proj_probs{marg_id} = obj.proj_probs{marg_id} / sum(obj.proj_probs{marg_id});
            end
        end

        function setRandStream(obj, rand_stream)
            % Setter of the rand_stream memeber variable

            obj.rand_stream = rand_stream;
        end
        
        function samps = randSample(obj, samp_num)
            % Randomly generate samples from the worst-case probability measure
            % Input:
            %   samp_num: number of samples
            % Output:
            %   samps: matrix containing the generated samples as rows

            samp_dual_point_indices = randsample(obj.rand_stream, length(obj.dual_points_probs), samp_num, true, ...
                obj.dual_points_probs);

            samps = zeros(samp_num, obj.problem.unc_num);

            for marg_id = 1:obj.problem.unc_num
                samp_proj_indices = obj.proj_mappings{marg_id}(samp_dual_point_indices);
                samps(:, marg_id) = obj.problem.uncertain.marginals{marg_id}.randSampleConditional( ...
                    obj.proj_atoms{marg_id}, obj.proj_probs{marg_id}, samp_proj_indices, obj.rand_stream);
            end
        end
    end
end

