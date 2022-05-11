load('exp/DRO_scheduling/exp_inputs.mat');

mixtrnorm_invcdf = cell(2, 1);

cdf_granularity = 5;

mu = mixtrnorm_mu;
sig = sqrt(mixtrnorm_sig2);
w = mixtrnorm_w;

cdf_bound1 = normcdf((mixtrnorm_trunc(1) - mu) ./ sig);
cdf_bound2 = normcdf((mixtrnorm_trunc(2) - mu) ./ sig);
normconst = cdf_bound2 - cdf_bound1;

t_func = @(v)(sum(((normcdf(((v' - mu) ./ sig)) - cdf_bound1) ...
    ./ normconst) .* w, 1)');

xx = (0:10^cdf_granularity)' / 10^cdf_granularity ...
    * (mixtrnorm_trunc(2) - mixtrnorm_trunc(1)) + mixtrnorm_trunc(1);
yy = t_func(xx);

zero_id = find(yy == 0, 1, 'last');
yy = yy(zero_id:end);
xx = xx(zero_id:end);

[yy, ia] = unique(yy);
xx = xx(ia);

mixtrnorm_invcdf{1} = yy;
mixtrnorm_invcdf{2} = xx;

save('exp/DRO_scheduling/exp_invcdf.mat', 'mixtrnorm_invcdf', ...
    '-v7.3');