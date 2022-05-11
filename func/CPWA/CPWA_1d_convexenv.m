function [conv_knots, conv_vals, conv_vals_old] ...
    = CPWA_1d_convexenv(knots, vals)
% Compute the convex envelope of a one-dimensional piece-wise affine
% function with compact domain over a compact interval
% Inputs: 
%       knots: the knots in the one-dimensional piece-wise affine function
%       (the first and last knots define the domain)
%       vals: the corresponding values of the one-dimensional piece-wise
%       affine function at the knots
% Outputs: 
%       conv_knots: the knots in the convex envelope function
%       conv_vals: the corresponding values of the convex envelope function
%       conv_vals_old: the values of the convex envelope function at the
%       input knots

% check the inputs
knot_no = length(knots);
assert(length(vals) == knot_no, 'input function mis-specified');

conv_knots = zeros(knot_no, 1);
conv_vals = zeros(knot_no, 1);
conv_vals_old = zeros(knot_no, 1);

% the first knot of the envelope is the same as the input function
conv_knots(1) = knots(1);
conv_vals(1) = vals(1);
conv_vals_old(1) = vals(1);

% the current knot touching the convex envelope function
cur_knot = 1;

% the counter for the number of knots in the convex envelope function
counter = 1;

% begin the loop
while cur_knot < knot_no
    % compute the slope of all the potential supporting hyperplanes
    slopes = (vals(cur_knot + 1:end) - vals(cur_knot)) ...
        ./ (knots(cur_knot + 1:end) - knots(cur_knot));
    
    % choose the lowest supporting hyperplane
    slope_min = min(slopes);
    
    % allow a small numerical inaccuracy to avoid having "almost flat"
    % segments in the convex envelop
    min_id = find(slopes <= slope_min + 1e-10, 1, 'last');

    % evaluate the convex envelope function at the intermediate knots
    conv_vals_old(cur_knot + 1:cur_knot + min_id - 1) = vals(cur_knot) ...
        + (vals(cur_knot + min_id) - vals(cur_knot)) ...
        * (knots(cur_knot + 1:cur_knot + min_id - 1) - knots(cur_knot)) ...
        / (knots(cur_knot + min_id) - knots(cur_knot));
    conv_vals_old(cur_knot + min_id) = vals(cur_knot + min_id);
    
    % add the new knot 
    cur_knot = cur_knot + min_id;
    counter = counter + 1;
    conv_knots(counter) = knots(cur_knot);
    conv_vals(counter) = vals(cur_knot);
end

conv_knots = conv_knots(1:counter);
conv_vals = conv_vals(1:counter);

end