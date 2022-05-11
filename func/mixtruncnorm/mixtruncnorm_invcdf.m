function [z, iter] = mixtruncnorm_invcdf(x, mu, sig2, w, bounds, tol)
% Compute the inverse cumulative distribution function (cdf) of the
% mixture of truncated normal distributions using bisection
% Inputs:
%       x: inputs to the inverse cdf
%       mu: mu parameter of each mixture component
%       sig2: sigma^2 parameter of each mixture component
%       w: weight of each mixture compoment
%       bounds: the lower and upper truncation points (same for all the
%       components)
%       tol: the numerical tolerance (default is 1e-8)
% Outputs:
%       z: the corresponding inverse cdf values
%       iter: number of iterations

n = length(x);

if ~exist('tol', 'var') || isempty(tol)
    tol = 1e-8;
end

sig = sqrt(sig2);

% the normalizing constants of each component
cdf_bound1 = normcdf((bounds(1) - mu) ./ sig);
cdf_bound2 = normcdf((bounds(2) - mu) ./ sig);
normconst = cdf_bound2 - cdf_bound1;

if length(w) == 1
    % if there is only a single component, then norminv can be used
    z = norminv(x * normconst + cdf_bound1) * sqrt(sig2) + mu;
    iter = 0;
else
    % the cdf
    t_func = @(v)(sum(((normcdf(((v' - mu) ./ sig)) - cdf_bound1) ...
        ./ normconst) .* w, 1)');

    % initial interval
    z_range = repmat(bounds', n, 1);
    f_range = [t_func(z_range(:, 1)), t_func(z_range(:, 2))];
    iter = 0;

    while max(f_range(:, 2) - f_range(:, 1)) > tol
        z_mid = mean(z_range, 2);
        f_mid = t_func(z_mid);
        left_list = f_mid <= x;
        right_list = f_mid > x;
        z_range(left_list, 1) = z_mid(left_list);
        z_range(right_list, 2) = z_mid(right_list);
        f_range = [t_func(z_range(:, 1)), t_func(z_range(:, 2))];
        iter = iter + 1;
    end

    z = mean(z_range, 2);
end

end