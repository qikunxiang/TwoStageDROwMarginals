function [bd, bd_list] = mixtruncnorm_wassrad_bounded( ...
    mu, sig2, w, bounds, knots)
% Compute a bound on the Wassertein-1 radius of the moment set defined by
% one-dimensional CPWA functions
% Inputs:
%       mu: mu parameter of each mixture component
%       sig2: sigma^2 parameter of each mixture component
%       w: weight of each mixture component
%       bounds: the end points of the support (same for all the components)
%       knots: the knots defining the alternative basis CPWA functions
% Outputs: 
%       bd: the bound on radius
%       bd_list: the list of bounds on the terms corresponding to each
%       interval

% check the inputs
compno = length(mu);

assert(length(sig2) == compno, 'mixture components mis-specified');
assert(length(w) == compno, 'mixture components mis-specified');
assert(all(knots > bounds(1) & knots < bounds(2)), ...
    ['all knots must be strictly between the lower and upper ', ...
    'truncation points']);

% add two additional knots correspoinding to the end points
knots = [bounds(1); knots; bounds(2)];

knots_diff = diff(knots);
assert(all(knots_diff > 0), 'all knots must be strictly increasing');

% compute the moment values
momts = mixtruncnorm_partialexp(mu, sig2, w, bounds, knots, bounds(2), ...
    1, -knots);

% compute upper and lower bounds on the cdf of the unknown measure
cdf_bounds = max(0, min(1, (knots_diff + diff(momts)) ...
    ./ knots_diff));
cdf_bounds(1) = 0;
cdf_bounds(end) = 1;

% the integral limits of the terms in the sum
int_lims = mixtruncnorm_invcdf(cdf_bounds, mu, sig2, w, bounds);

% the mid-points of the (j-1)-th and (j+1)-th knots
mid_pts = (knots(1:end - 2) + knots(3:end)) / 2;

% the first half of the integral
v1 = mixtruncnorm_partialexp(mu, sig2, w, bounds, int_lims(1:end - 1), ...
    min(mid_pts, int_lims(2:end)), -1, knots(3:end));

% the second half of the integral
v2 = mixtruncnorm_partialexp(mu, sig2, w, bounds, ...
    max(mid_pts, int_lims(1:end - 1)), ...
    int_lims(2:end), 1, -knots(1:end - 2));

bd_list = v1 + v2;

bd = sum(bd_list);

end