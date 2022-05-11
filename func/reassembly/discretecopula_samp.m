function [samp, samp_atom_id] = discretecopula_samp(atoms, probs, ...
    samp_no, rand_stream, samp_mode)
% Generate random samples from the copula of a given discrete distribution
% Inputs: 
%       atoms: atoms in the discrete joint distribution, each row
%       represents an atom
%       probs: the corresponding probabilities at the atoms
%       samp_no: number of samples
%       rand_stream: RandStream object to guarantee reproducibility
%       samp_mode: the copula used for generating uniform random variables,
%       either 'comonotone' or 'independent' (default is 'comonotone')
% Outputs: 
%       samp: samples, each row represents a sample
%       samp_atom_id: samples of the associated atom indices

if ~exist('samp_mode', 'var') || isempty(samp_mode)
    samp_mode = 'comonotone';
end

[atom_no, d] = size(atoms);

[~, temp] = sort(atoms);
probs_marg = probs(temp);
[~, atoms_order] = sort(temp);
atoms_order_lin = atoms_order + (0:d - 1) * atom_no;
probs_cum = cumsum(probs_marg, 1);
atoms_unit = probs_cum(atoms_order_lin);
probs_rem = probs_marg(atoms_order_lin);

samp_atom_id = randsample(rand_stream, atom_no, samp_no, true, probs);

if strcmp(samp_mode, 'comonotone')
    unif_sample = repmat(rand(rand_stream, samp_no, 1), 1, d);
elseif strcmp(samp_mode, 'independent')
    unif_sample = rand(rand_stream, samp_no, d);
else
    error('unknown sampling mode');
end

samp = atoms_unit(samp_atom_id, :) - unif_sample ...
    .* probs_rem(samp_atom_id, :);

samp = max(min(samp, 1), 0);

end

