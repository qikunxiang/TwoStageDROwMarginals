function v = norm_partialexp(mu, sig2, k1, k2, a, b)
% Compute the partial expectation of a normal distribution of the
% form E[(aX+b)I{k1,k2}], where mu, sig2, a, b, k1, k2 are vectors of the
% same length (or constants)
% Inputs:
%       mu: the mean parameter of the normal distribution
%       sig2: the variance parameter (sigma squared) of the normal
%           distribution
%       k1: the left end point of the interval (lower integration limit)
%       k2: the right end point of the interval (upper integration limit)
%       a: the coefficient of X
%       b: the intercept
% Outputs:
%       v: the value of the partial expectation

if length(k1) == 1 && length(k2) > 1
    k1 = repmat(k1, length(k2), 1);
elseif length(k2) == 1 && length(k1) > 1
    k2 = repmat(k2, length(k1), 1);
end

k2 = max(k1, k2);

% standardize the end points
mat = ([k1, k2] - mu) ./ sqrt(sig2);

P = normcdf(mat);
expmat = exp(-0.5 * (mat .^ 2));

v = (a .* mu + b) .* (P(:, 2) - P(:, 1)) + (1 / sqrt(2 * pi)) ...
    * a .* sqrt(sig2) .* (expmat(:, 1) - expmat(:, 2));

end