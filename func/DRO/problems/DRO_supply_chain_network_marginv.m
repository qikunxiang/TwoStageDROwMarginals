function samps = DRO_supply_chain_network_marginv(i, u, ...
    customer_no, invcdf_cell, failure_prob)
% Utility function for computing the inverse cdf transform in the supply
% chain network design problem with edge failure
% Inputs: 
%       i: the dimension index
%       u: vector of (non-independent) uniformly distributed random
%       variables
%       customer_no: number of customers
%       invcdf_cell: cell array containing the information about the
%       inverse cdf of the demand distributions
%       failure_prob: vector containing the failure probability of the
%       edges
% Output:
%       samps: matrix containing generated samples in rows

if i <= customer_no
    samps = interp1(invcdf_cell{i, 1}, invcdf_cell{i, 2}, u);
else
    samps = u > (1 - failure_prob(i - customer_no));
end

end

