function v = mixtruncnorm_partialexp(mu, sig2, w, bounds, k1, k2, a, b)
% Compute the partial expectation of a mixture of truncated normal
% distributions of the form E[(aX+b)I{k1,k2}], where mu, sig2, w are
% vectors of the same length specifying the mean, variance, and weight of
% each component of the mixture, bounds specify the truncation points, k1,
% k2, a, b are vectors of the same length (or scalars)
% Inputs:
%       mu: the mean parameter of each component
%       sig2: the variance parameter (sigma squared) of each component
%       w: the weight of each component
%       bounds: the lower and upper truncation points (same for all
%       components)
%       k1: the left end point of the interval (lower integration limit)
%       k2: the right end point of the interval (upper integration limit)
%       a: the coefficient of X
%       b: the intercept
% Outputs:
%       v: the value of the partial expectation

sig = sqrt(sig2);
t1 = bounds(1);
t2 = bounds(2);

% the normalizing constants of each component
normconst = normcdf((t2 - mu) ./ sig) - normcdf((t1 - mu) ./ sig);

% truncate the integration limits if necessary
k1 = max(k1, t1);
k2 = min(k2, t2);

v = [];

% add up the contribution from all components
for comp = 1:length(mu)
    vcomp = norm_partialexp(mu(comp), sig2(comp), k1, k2, a, b) ...
        / normconst(comp) * w(comp);
    
    if isempty(v)
        v = vcomp;
    else
        v = v + vcomp;
    end
end

end