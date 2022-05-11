function [atoms, prob] = momentset1d_measinit_bounded(knots, expval)
% Compute a discrete measure that satisfies all moment constraints with
% CPWA functions
% Inputs: 
%       knots: the knots defining the alternative basis CPWA functions
%       where the first and last knots represent the end points of the
%       closed interval
%       expval: the expected value of the alternative basis CPWA functions
% Outputs: 
%       atoms: the position of the atoms
%       prob: the corresponding probabilities of the atoms


atoms = knots;
prob = expval;

end

