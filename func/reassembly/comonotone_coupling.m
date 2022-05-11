function [atoms, probs] = comonotone_coupling(marg)
% Compute the comonotone coupling of finitely-supported marginal
% distributions
% Inputs: 
%       marg: N-by-2 cell array where the first column contains the atoms
%       of the marginal distribution and the second column contains the
%       corresponding probabilities
% Outputs: 
%       atoms: the list of atoms in the joint distribution
%       probs: the corresponding probabilities of the atoms

N = size(marg, 1);

atomno_list = zeros(N, 1);

for i = 1:N
    % check the inputs
    assert(length(marg{i, 1}) == length(marg{i, 2}), ...
        'input mis-specified');
    
    % normalize the probabilities (in case there is a small numerical
    % error) 
    marg{i, 2} = marg{i, 2} / sum(marg{i, 2});
    
    nonzero_list = marg{i, 2} > 0;
    marg{i, 2} = marg{i, 2}(nonzero_list);
    marg{i, 1} = marg{i, 1}(nonzero_list);
    
    atomno_list(i) = length(marg{i, 1});
end

atomno_max = max(atomno_list);

% place atoms and probabilities in columns of matrices
atom_mat = zeros(atomno_max, N);
prob_mat = zeros(atomno_max, N);

for i = 1:N
    atom_mat(1:atomno_list(i), i) = marg{i, 1};
    prob_mat(1:atomno_list(i), i) = marg{i, 2};
end

% initialization before the loop
joint_atomno = 0;
joint_atom = zeros(nnz(prob_mat), N);
joint_prob = zeros(nnz(prob_mat), 1);
counter = sub2ind([atomno_max, N], ones(N, 1), (1:N)');
counter_max = sub2ind([atomno_max, N], atomno_list, (1:N)');

while true
    % take out one atom
    prob_min = min(prob_mat(counter));
    
    % record the atom and its probability
    joint_atomno = joint_atomno + 1;
    joint_atom(joint_atomno, :) = atom_mat(counter);
    joint_prob(joint_atomno) = prob_min;
    
    % decrease the probability from the remaining probabilities
    prob_mat(counter) = prob_mat(counter) - prob_min;
    
    % advance the counter for the corresponding column
    counter = counter + (prob_mat(counter) == 0);
    
    if any(counter > counter_max)
        break;
    end
end

% prepare output
atoms = joint_atom(1:joint_atomno, :);
probs = joint_prob(1:joint_atomno);

end