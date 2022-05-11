function [v, wrad_ub, int_probs] = mixtruncnorm_momentset(mu, sig2, w, ...
    knots)
% Construct a moment set around a mixture of truncated normal distribution
% with bounded support
% Inputs:
%       mu: mu parameter of each mixture component
%       sig2: sigma^2 parameter of each mixture component
%       w: weight of each mixture compoment
%       knots: knots of the CPWA functions where the first and last knots
%       represent the points of truncation (same for all the components)
% Output: 
%       v: values of the expectation of the basis functions
%       wrad_ub: a (crude) upper bound for the W1 radius of the resulting
%       moment set
%       int_probs: the probability in each interval between two knots

t1 = knots(1);
t2 = knots(end);

% check the inputs
compno = length(mu);

assert(length(sig2) == compno, 'mixture components mis-specified');
assert(length(w) == compno, 'mixture components mis-specified');

k1 = knots(1:end - 1);
k2 = knots(2:end);
kdiff = k2 - k1;
assert(all(kdiff > 0), 'knots must be in ascending order');

% integral to the right of each knot
pexp1 = mixtruncnorm_partialexp(mu, sig2, w, [t1; t2], ...
    k1, k2, -1 ./ kdiff, k2 ./ kdiff);
% integral to the left of each knot
pexp2 = mixtruncnorm_partialexp(mu, sig2, w, [t1; t2], ...
    k1, k2, 1 ./ kdiff, -k1 ./ kdiff);

% compute the probability in each interval
int_probs = pexp1 + pexp2;

% add up the integrals
v = [pexp1(1); pexp1(2:end) + pexp2(1:end - 1); pexp2(end)];

% compute an upper bound for the W1 radius
wrad_ub = 2 * sum(int_probs .* kdiff);

end