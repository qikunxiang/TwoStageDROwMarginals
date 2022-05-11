load('exp/DRO_assembly/exp_inputs.mat');

mixtrnorm_invcdf_cell = cell(N, 2);
cdf_granularity = 5;

for i = 1:N
    mu = mixtrnorm_mu_cell{i};
    sig = sqrt(mixtrnorm_sig2_cell{i});
    w = mixtrnorm_w_cell{i};

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

    mixtrnorm_invcdf_cell{i, 1} = yy;
    mixtrnorm_invcdf_cell{i, 2} = xx;
end

save('exp/DRO_assembly/exp_invcdf.mat', 'mixtrnorm_invcdf_cell', ...
    '-v7.3');