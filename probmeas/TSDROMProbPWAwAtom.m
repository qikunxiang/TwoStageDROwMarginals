classdef TSDROMProbPWAwAtom < TSDROMProb
    % Class for probability measures supported on a one-dimensional interval which is a mixture of a discrete measure and a continuous
    % measure with (possibly discontinuous) piece-wise affine (PWA) density
    
    properties(GetAccess = public, SetAccess = protected)
        Dens = struct;
    end
    
    methods(Access = public)
        function obj = TSDROMProbPWAwAtom(pwa_intv, pwa_dens, disc_atoms, disc_probs)
            % Constructor
            % The lower and upper bound of the support interval will be the smallest interval containing all the atoms and
            % subintervals. The CPA part of the measure will be renormalized to (1 - sum(disc_probs)). The subintervals in the PWA part
            % must be non-overlapping.
            %   Inputs:
            %       pwa_intv: two-column matrix containing the left and right end points of the subintervals in the PWA part
            %       pwa_dens: two-column matrix containing the densities at the left and right end points of the subintervals in the 
            %       PWA part
            %       disc_atoms: vector containing the atoms in the discrete part
            %       disc_probs: vector containing the probabilities of the atoms in the discrete part
            
            % check that the subintervals are non-empty
            assert(all(pwa_intv(:, 2) > pwa_intv(:, 1)), 'some subintervals are empty');
            assert(size(pwa_intv, 2) == 2, 'subintervals mis-specified');
            intv_num = size(pwa_intv, 1);

            % check that the subintervals are non-overlapping
            [~, intv_sorted_indices] = sort(pwa_intv(:, 1), 'ascend');
            pwa_intv = pwa_intv(intv_sorted_indices, :);
            pwa_dens = pwa_dens(intv_sorted_indices, :);
            assert(all(pwa_intv(1:end - 1, 2) <= pwa_intv(2:end, 1)), 'some subintervals are overlapping');

            % check that the densities are correctly specified
            assert(size(pwa_dens, 1) == intv_num && size(pwa_dens, 2) == 2, 'subinterval densities mis-specified');
            assert(all(all(pwa_dens >= 0)), 'some densities are negative');

            % check that the atom probabilities are correctly specified
            assert(isvector(disc_atoms), 'discrete atoms mis-specified');
            disc_atom_num = length(disc_atoms);
            assert(size(disc_probs, 1) == disc_atom_num && size(disc_probs, 2) == 1, 'discrete atom probabilities mis-specified');

            % check that the atoms have positive probabilities
            assert(all(disc_probs >= 1e-8), 'some discrete atoms have negligible probabilities');

            [disc_atoms, disc_sorted_indices] = sort(disc_atoms, 'ascend');
            assert(all(diff(disc_atoms) > 1e-8), 'some discrete atoms are almost non-distinct');
            disc_probs = disc_probs(disc_sorted_indices);

            assert(intv_num > 0 || disc_atom_num > 0, 'the probability measure is empty');

            % set the lower and upper bounds of the support interval
            if intv_num > 0 && disc_atom_num > 0
                intv_lb = min(pwa_intv(1, 1), disc_atoms(1));
                intv_ub = max(pwa_intv(end, 2), disc_atoms(end));
            elseif intv_num > 0
                intv_lb = pwa_intv(1, 1);
                intv_ub = pwa_intv(end, 2);
            else
                intv_lb = disc_atoms(1);
                intv_ub = disc_atoms(end);
            end

            obj@TSDROMProb(intv_lb, intv_ub);

            % normalize the PWA part of the density
            disc_total_probs = sum(disc_probs);
            pwa_total_probs = sum(sum(pwa_dens, 2) .* (pwa_intv(:, 2) - pwa_intv(:, 1)) / 2);
            pwa_dens = pwa_dens .* ((1 - disc_total_probs) / pwa_total_probs);
            
            assert(~any(all(pwa_dens < 1e-8, 2)), 'some subintervals have negligible probabilities');

            % compute the slope and intercept of the affine functions in the subintervals
            pwa_dens_slp = (pwa_dens(:, 2) - pwa_dens(:, 1)) ./ (pwa_intv(:, 2) - pwa_intv(:, 1));
            pwa_dens_itc = pwa_dens(:, 1) - pwa_intv(:, 1) .* pwa_dens_slp;

            % set the member variables
            obj.Dens.pwa_intv = pwa_intv;
            obj.Dens.pwa_dens = pwa_dens;
            obj.Dens.pwa_intv_num = intv_num;
            obj.Dens.pwa_dens_slp = pwa_dens_slp;
            obj.Dens.pwa_dens_itc = pwa_dens_itc;
            obj.Dens.disc_atoms = disc_atoms;
            obj.Dens.disc_probs = disc_probs;
            obj.Dens.disc_atom_num = disc_atom_num;

            obj.generateInternalCDFRep();
        end

        function x = computeLContInvCDF(obj, p)
            % Compute the left-continuous generalized inverse cumulative density function
            % Input:
            %   p: vector containing the probabilities with respect to which the inverse CDF is to be evaluated
            % Output:
            %   x: vector containing the computed inverse CDF

            p = max(min(p, 1), 0);

            % select the subinterval or atom in the internal CDF representation that each probability value belongs to
            rep = obj.Dens.cdf_rep;
            [~, elem_indices] = max(p - rep(:, 6)' <= 3e-16, [], 2);

            % solve a quadratic equation in the corresponding subinterval or atom to determine the inverted value
            rep = rep(elem_indices, :);
            prob_ratio = (p - rep(:, 6)) ./ rep(:, 5) + 1;
            prob_left = rep(:, 3);
            prob_right = rep(:, 4);
            x_ratio = prob_ratio .* (prob_left + prob_right) ./ (sqrt(prob_ratio .* prob_right.^2 ...
                + (1 - prob_ratio) .* prob_left.^2) + prob_left);
            x = rep(:, 1) + x_ratio .* (rep(:, 2) - rep(:, 1));
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

            % integrals with respect to the piece-wise affine part
            int_lb = max(lb, obj.Dens.pwa_intv(:, 1)');
            int_ub = max(min(ub, obj.Dens.pwa_intv(:, 2)'), int_lb);
            int_cube_coef = a .* obj.Dens.pwa_dens_slp' / 3;
            int_quad_coef = (a .* obj.Dens.pwa_dens_itc' + b .* obj.Dens.pwa_dens_slp') / 2;
            int_line_coef = b .* obj.Dens.pwa_dens_itc';

            int_intv = int_cube_coef .* (int_ub.^3 - int_lb.^3) + int_quad_coef .* (int_ub.^2 - int_lb.^2) ...
                + int_line_coef .* (int_ub - int_lb);

            % integrals with respect to the discrete part
            int_disc = (a .* obj.Dens.disc_atoms' + b) .* (lb <= obj.Dens.disc_atoms' & ub > obj.Dens.disc_atoms') ...
                .* obj.Dens.disc_probs';

            vals = sum(int_intv, 2) + sum(int_disc, 2);
        end
    end

    methods(Access = protected)
        function generateInternalCDFRep(obj)
            % Generate an internal representation of the cumulative distribution function, in which the subintervals and atoms are
            % arranged from left to right

            % the representation is a matrix with 6 columns containing both the subintervals with piece-wise affine (PWA) densities and
            % the atoms;
            % a subinterval is split if atoms are present within it;
            % column 1 contains the left end points of the subintervals;
            % column 2 contains the right end points of the subintervals; atoms correspond to subintervals where the two end points are
            % identical;
            % column 3 contains the densities at the left end points;
            % column 4 contains the densities at the right end points; for an atom, the densities on both the left and right end points
            % are equal to the atom probability;
            % column 5 contains the total probabilities in the subintervals;
            % column 6 contains the cumulative probabilities to the left of the subintervals

            rep = zeros(obj.Dens.pwa_intv_num + 2 * obj.Dens.disc_atom_num, 6);
            atom_indices = false(obj.Dens.pwa_intv_num + 2 * obj.Dens.disc_atom_num, 1);
            counter = 0;
            
            % create a matrix with 3 columns: column 1 contains the left end points of the subintervals and the atoms, column 2
            % contains the subinterval/atom indices; column 3 contains logical indicators for whether a row represents a subinterval or
            % an atom
            union = [obj.Dens.pwa_intv(:, 1), (1:obj.Dens.pwa_intv_num)', false(obj.Dens.pwa_intv_num, 1); ...
                obj.Dens.disc_atoms, (1:obj.Dens.disc_atom_num)', true(obj.Dens.disc_atom_num, 1)];

            [~, sorted_indices] = sort(union(:, 1), 'ascend');
            union = union(sorted_indices, :);

            % scan through the rows of this matrix and create subinterval representations in the rep matrix
            for row_id = 1:size(union, 1)
                left_end = union(row_id, 1);
                elem_id = union(row_id, 2);
                is_atom = union(row_id, 3);

                counter = counter + 1;
                atom_indices(counter) = is_atom;
                rep(counter, 1) = left_end;

                if is_atom
                    rep(counter, 2) = left_end;
                    rep(counter, 3:4) = obj.Dens.disc_probs(elem_id);
                else
                    rep(counter, 2) = obj.Dens.pwa_intv(elem_id, 2);
                    rep(counter, 3:4) = obj.Dens.pwa_dens(elem_id, :);
                end

                if ~is_atom || counter == 1 || atom_indices(counter - 1)
                    continue;
                end

                % the only case that would require special treatment is when an atom is added and the previous element in the
                % representation is a subinterval
                if left_end < rep(counter - 1, 1) + 1e-8
                    % if the current atom is very close to the left end point of the previous subinterval, swap the atom with the
                    % previous subinterval
                    rep(counter - 1, 1) = left_end;
                    temp_row = rep(counter - 1, :);
                    rep(counter - 1, :) = rep(counter, :);
                    rep(counter, :) = temp_row;
                    atom_indices(counter) = false;
                    atom_indices(counter - 1) = true;
                elseif left_end < rep(counter - 1, 2) - 1e-8
                    % if the current atom is contained in the previous subinterval, split the previous subinterval into two halfs
                    ratio = (left_end - rep(counter - 1, 1)) / (rep(counter - 1, 2) - rep(counter - 1, 1));
                    dens_mid = rep(counter - 1, 3) + ratio * (rep(counter - 1, 4) - rep(counter - 1, 3));
                    rep(counter + 1, 1) = left_end;
                    rep(counter + 1, 2) = rep(counter - 1, 2);
                    rep(counter + 1, 3) = dens_mid;
                    rep(counter + 1, 4) = rep(counter - 1, 4);
                    rep(counter - 1, 2) = left_end;
                    rep(counter - 1, 4) = dens_mid;
                    atom_indices(counter + 1) = false;
                    counter = counter + 1;
                end
            end

            rep = rep(1:counter, :);
            atom_indices = atom_indices(1:counter);

            rep(:, 5) = (rep(:, 3) + rep(:, 4)) .* (rep(:, 2) - rep(:, 1)) / 2;
            rep(atom_indices, 5) = rep(atom_indices, 3);
            rep(:, 6) = cumsum(rep(:, 5));

            % fix the cumulative probability of the last element to be 1 to prevent numerical issues
            rep(end, 6) = 1;

            obj.Dens.cdf_rep = rep;
        end
    end
end

