function [A, rhs] = DRO_feascons(x, lambda, knot_vec, knot_offset, ...
    keep_list, DRO)
% Generate linear feasibility constraints in the LSIP problem for solving a
% two-stage distributionally robust optimization (DRO) problem with
% marginal constraints
% Inputs: 
%       x: each row is a set of knot indices
%       lambda: each row is a set of values of the second-stage decision
%       knot_vec: the knots in each dimension placed into a column vector
%       knot_offset: the offset of each dimension in the column vector
%       keep_list: list of indices to keep intended for removing the first
%       know from dimensions 2 to N for identification
%       DRO: the specification of the DRO problem
% Outputs: 
%       A, rhs: the resulting linear constraints A * [a; y] >= rhs

constr_no = size(x, 1);
N = size(x, 2);

x_w_offset = x + knot_offset';
x_val = knot_vec(x_w_offset);

if constr_no == 1
    x_val = reshape(x_val, 1, []);
end

rhs = full(sum((x_val * DRO.stage2.W' + DRO.stage2.b') .* lambda, 2));

% sanitize the output to prevent potential numerical issues
rhs(abs(rhs) < 1e-9) = 0;

A_y = sparse(repmat((1:constr_no)', N, 1), x_w_offset(:), 1, ...
    constr_no, length(knot_vec));

A = [-sparse(lambda * DRO.stage2.V), A_y(:, keep_list)];

% sanitize the output to prevent potential numerical issues
A(abs(A) < 1e-9) = 0;

end

