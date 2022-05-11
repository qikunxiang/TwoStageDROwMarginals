function [atoms_new, probs_new] = reassembly_increment(marg_atoms_new, ...
    marg_probs_new, atoms, probs)
% Compute the reassembly of a discrete measure with discrete 1D marginals
% Inputs: 
%       marg_atoms_new: cell array containing atoms in each marginal
%       marg_probs_new: cell array containing probabilities in of atoms in
%       each marginal
%       atoms: the atoms in the input measure
%       probs: the corresponding probabilities in the input measure
% Outputs: 
%       atoms_new: the atoms in the reassembly
%       probs_new: the corresponding probabilities in the reassembly

N = length(marg_atoms_new);

% check the inputs
assert(size(atoms, 2) == N, 'atoms mis-specified');
assert(size(atoms, 1) == length(probs), 'probabilities mis-specified');

% normalize the probabilities (in case there is a small numerical error)
probs = probs / sum(probs);

for i = 1:N
    assert(length(marg_probs_new{i}) == length(marg_atoms_new{i}), ...
        'new probabilities mis-specified');
end

% compute the updated measure

% remove atoms with zero probability
nonzero_list = probs > 0;
probs_new = probs(nonzero_list);
atoms_new = atoms(nonzero_list, :);

% begin the loop to compute a reassembly
for i = 1:N
    marg_atoms_i = marg_atoms_new{i};
    marg_probs_i = marg_probs_new{i};
    
    % sort the atoms into ascending order in the i-th dimension
    [~, asc_order] = sort(atoms_new(:, i), 'ascend');
    atoms_sorted = atoms_new(asc_order, :);
    prob_sorted = probs_new(asc_order);
    
    % the atoms and probabilities in the glued distribution
    atomno_max = size(atoms_sorted, 1) * length(marg_atoms_i);
    atoms_glue = zeros(atomno_max, N + 1);
    prob_glue = zeros(atomno_max, 1);
    
    % initialize the counters
    counter_joint = 1;
    counter_marg = 1;
    counter_glue = 0;
    
    % begin the loop to couple the marginals
    while true
        % take out one atom
        prob_min = min(prob_sorted(counter_joint), ...
            marg_probs_i(counter_marg));
        
        % record the atom and its probability
        counter_glue = counter_glue + 1;
        atoms_glue(counter_glue, :) = [atoms_sorted(counter_joint, :), ...
            marg_atoms_i(counter_marg)];
        prob_glue(counter_glue) = prob_min;
        
        % decrease the probability from the remaining probabilities
        prob_sorted(counter_joint) = prob_sorted(counter_joint) ...
            - prob_min;
        marg_probs_i(counter_marg) = marg_probs_i(counter_marg) ...
            - prob_min;
        
        % advance the counters
        if prob_sorted(counter_joint) == 0
            counter_joint = counter_joint + 1;
        end
        
        if marg_probs_i(counter_marg) == 0
            counter_marg = counter_marg + 1;
        end
        
        if counter_joint > size(atoms_sorted, 1) ...
                || counter_marg > length(marg_atoms_i)
            break;
        end
    end
    
    % update the distribution
    atoms_glue = atoms_glue(1:counter_glue, :);
    atoms_new = atoms_glue(:, 1:N);
    atoms_new(:, i) = atoms_glue(:, end);
    probs_new = prob_glue(1:counter_glue);
end

end

