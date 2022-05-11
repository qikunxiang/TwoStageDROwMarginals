function [knots, knots_hist, bd_hist] ...
    = mixtruncnorm_momentset_greedy(mu, sig2, w, bounds, n)
% Greedily search for a moment set around a mixture of truncated normal
% distributions defined by one-dimensional CPWA functions with a given
% number of knots
% Inputs: 
%       mu: mu parameter of each mixture component
%       sig2: sigma^2 parameter of each mixture component
%       w: weight of each mixture component
%       bounds: the end points of the support (same for all the components)
%       n: number of knots required (at least 3)
% Outputs: 
%       knots: the resulted knots
%       knots_hist: the resulted knots in the order they are added
%       bd_hist: the history of the bounds

assert(n >= 3, 'n is too small');

% start with one knot at the quartiles
knots_int = mixtruncnorm_invcdf(0.5, mu, sig2, w, bounds);

% the knots and the boundary points
knots = [bounds(1); knots_int; bounds(2)];

% the values of the cdf at the knots
cdf_vals = [0; 0.5; 1];

knots_hist = [knots; zeros(n - length(knots), 1)];
bd_hist = inf(n, 1);
[~, bd_hist(length(knots))] = mixtruncnorm_momentset(mu, sig2, w, knots);

while length(knots) < n
    best_added_knot_cdf = (cdf_vals(1) + cdf_vals(2)) / 2;
    best_added_knot = mixtruncnorm_invcdf(best_added_knot_cdf, ...
        mu, sig2, w, bounds);
    best_cdf_vals = [cdf_vals(1); best_added_knot_cdf; cdf_vals(2:end)];
    best_knots = [knots(1); best_added_knot; knots(2:end)];
    [~, best_bd] = mixtruncnorm_momentset(mu, sig2, w, best_knots);

    for j = 2:length(knots) - 1
        added_knot_cdf = (cdf_vals(j) + cdf_vals(j + 1)) / 2;
        added_knot = mixtruncnorm_invcdf(added_knot_cdf, ...
            mu, sig2, w, bounds);

        added_cdf_vals = [cdf_vals(1:j); added_knot_cdf; ...
            cdf_vals(j + 1:end)];
        trial_knots = [knots(1:j); added_knot; knots(j + 1:end)];
        [~, trial_bd] = mixtruncnorm_momentset(mu, sig2, w, trial_knots);

        if trial_bd < best_bd
            best_added_knot = added_knot;
            best_cdf_vals = added_cdf_vals;
            best_knots = trial_knots;
            best_bd = trial_bd;
        end
    end
    
    cdf_vals = best_cdf_vals;
    knots = best_knots;
    knots_hist(length(knots)) = best_added_knot;
    bd_hist(length(knots)) = best_bd;
end

end