function v = momentset1d_basisfunc_bounded(x, knots)
% Basis functions defining moment sets in 1D when the domain is a bounded
% closed interval
% Inputs: 
%       x: the inputs to the basis functions
%       knots: knots of the CPWA functions where the first and last knots
%       represent the end points of the closed interval
% Outputs: 
%       v: values of the interpolation basis functions evaluated at x

assert(all(x >= knots(1) & x <= knots(end)), 'inputs out of bound');

knots_diff = diff(knots);
assert(all(knots_diff > 0), 'all knots must be strictly increasing');


v_raw = max(x - knots', 0);

% compute the values of the basis functions
v_diff = -diff(v_raw, 1, 2);
v_diff_norm = v_diff ./ knots_diff';
v = [-diff(v_diff_norm, 1, 2), v_diff_norm(:, end)];

end