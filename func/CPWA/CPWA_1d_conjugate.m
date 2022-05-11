function [conj_knots, conj_vals] = CPWA_1d_conjugate(knots, vals, ...
    conj_range)
% Compute the convex conjugate of a one-dimensional convex piece-wise
% affine function with compact domain over a compact interval
% Inputs: 
%       knots: the knots in the one-dimensional convex piece-wise affine
%       function (the first and last knots define the domain)
%       vals: the corresponding values of the one-dimensional convex
%       piece-wise affine function at the knots; the convex function must
%       not contain any redundant knots
%       conj_range: the range over which the conjugate function is computed
%       (the domain of the conjugate function is the whole real line
%       because the input function has compact domain)
% Outputs: 
%       conj_knots: the knots in the convex conjugate function
%       conj_vals: the corresponding values of the convex conjugate
%       function

% check the inputs
knotno = length(knots);
assert(length(vals) == knotno, 'input function mis-specified');
slopes = diff(vals) ./ diff(knots);
assert(all(diff(slopes) >= 0), 'the input function is not convex');
assert(all(diff(slopes) > 0), ...
    'the input function contains redundant knots');

% the slope between knots in the input function corresponds to the knots in
% the conjugate function
conj_knots = slopes;
conj_vals = conj_knots .* knots(1:end - 1) - vals(1:end - 1);

if conj_range(1) >= conj_knots(end)
    % if the entire interval is to the right of the right-most knot
    conj_knots = conj_range;
    conj_vals = conj_range * knots(end) - vals(end);
elseif conj_range(2) <= conj_knots(1)
    % if the entire interval is to the left of the left-most knot
    conj_knots = conj_range;
    conj_vals = conj_range * knots(1) - vals(1);
else
    left_id = find(conj_knots > conj_range(1), 1, 'first');
    right_id = find(conj_knots < conj_range(2), 1, 'last');
    conj_knots = [conj_range(1); conj_knots(left_id:right_id); ...
        conj_range(2)];
    conj_vals = [conj_range(1) * knots(left_id) - vals(left_id); ...
        conj_vals(left_id:right_id); ...
        conj_range(2) * knots(right_id + 1) - vals(right_id + 1)];
end

end

