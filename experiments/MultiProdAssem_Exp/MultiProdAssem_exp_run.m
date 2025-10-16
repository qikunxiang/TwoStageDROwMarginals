CONFIG = MultiProdAssem_exp_config();

load(CONFIG.SAVEPATH_INPUTS);

marginals = cell(product_num, 1);

for marg_id = 1:product_num
    marginals{marg_id} = eval(sprintf('%s(marginals_info{marg_id, 2}{:})', marginals_info{marg_id, 1}));
end

uncertain = TSDROMUncertain(marginals);
problem = TSDROMMultiProdAssem(uncertain, part_costs, assembly_parts, product_revenues, part_salvage_prices);

LSIP_options = struct;
LSIP_options.dual_feasibility_tolerance = 1e-5;
LSIP_options.primal_prob_thres = 1e-7;
LSIP_options.OT_contrib_thres = 1e-6;
LSIP_options.dual_weight_thres = 1e-6;

LSIP_options.log_file = CONFIG.LOGPATH_LSIP_MAIN;
LSIP_options.display = true;
LSIP_options.reduce = struct;
LSIP_options.reduce.preserve_init_constr = true;
LSIP_options.reduce.min_slack = 1e-7;
LSIP_options.reduce.thres = 2e-3;
LSIP_options.reduce.thres_quantile = 0.8;
LSIP_options.reduce.max_iter = 1000;
LSIP_options.reduce.freq = 100;

LSIP_global_options = struct;
LSIP_global_options.display = true;
LSIP_global_options.log_file = CONFIG.LOGPATH_LSIP_GLOBAL;

LSIP_LP_options = struct;
LSIP_LP_options.FeasibilityTol = 1e-9;
LSIP_LP_options.Method = 2;
LSIP_LP_options.Presolve = 0;
LSIP_LP_options.Crossover = 0;
LSIP_LP_options.OutputFlag = 1;
LSIP_LP_options.LogToConsole = 1;
LSIP_LP_options.LogFile = CONFIG.LOGPATH_LSIP_LP;

NewVar_options = struct;
NewVar_options.display = true;
NewVar_options.log_file = CONFIG.LOGPATH_NEWVAR_MAIN;
NewVar_options.init_extrapolation = 40;
NewVar_options.max_extrapolation = 640;

NewVar_global_options = struct;
NewVar_global_options.PoolSolutions = 200;
NewVar_global_options.PoolSearchMode = 2;
NewVar_global_options.PoolGap = 0.95;
NewVar_global_options.ConcurrentMIP = 1;
NewVar_global_options.NodefileStart = 0.5;
NewVar_global_options.NodefileDir = '../GurobiNodes';
NewVar_global_options.OutputFlag = 1;
NewVar_global_options.LogToConsole = 1;
NewVar_global_options.LogFile = CONFIG.LOGPATH_NEWVAR_GLOBAL;

LSIP_solver = TSDROMDualLSIPSolver(LSIP_options, LSIP_LP_options, LSIP_global_options);
NewVar_solver = TSDROMNewVariableINCSolver(NewVar_options, NewVar_global_options);


Iter_options = struct;
Iter_options.display = true;
Iter_options.log_file = CONFIG.LOGPATH_ITERSOLVER;

Iter_options.phase1 = struct;
Iter_options.phase1.Kpot_tolerance = 100;
Iter_options.phase1.MIPGap_tolerance = 0.10;

Iter_options.phase2 = struct;

Iter_options.phase2.DualLSIPLPOptions = struct;
Iter_options.phase2.DualLSIPLPOptions.Method = -1;
Iter_options.phase2.DualLSIPLPOptions.Presolve = -1;
Iter_options.phase2.DualLSIPLPOptions.Crossover = -1;

Iter_options.phase2.NewVarGlobalOptions = struct;
Iter_options.phase2.NewVarGlobalOptions.PoolSolutions = 40;
Iter_options.phase2.NewVarGlobalOptions.PoolSearchMode = 0;

Iter_options.phase2.Kpot_tolerance_init = 100;
Iter_options.phase2.Kpot_tolerance_final = 10;
Iter_options.phase2.MIPGap_tolerance_init = 0.10;
Iter_options.phase2.MIPGap_tolerance_final = 0.05;
Iter_options.phase2.dual_points_redundant_threshold = 50;
Iter_options.phase2.iter_num = 1e4;

Iter_options.phase3 = struct;

Iter_options.phase3.NewVarGlobalOptions = struct;
Iter_options.phase3.NewVarGlobalOptions.TimeLimit = 300;

Iter_options.phase3.Kpot_tolerance = 1;
Iter_options.phase3.MIPGap_tolerance = 0.01;
Iter_options.phase3.iter_num = 200;

Iter_options.phase4 = struct;

Iter_options.phase4.NewVarGlobalOptions = struct;
Iter_options.phase4.NewVarGlobalOptions.TimeLimit = Inf;

Iter_options.phase4.Kpot_tolerance = 1;
Iter_options.phase4.MIPGap_tolerance = 0.01;

IterSolver = TSDROMIterativeSolver(LSIP_solver, NewVar_solver, Iter_options);
IterSolver.setProblem(problem);

tolerance = 10;
OT_tolerance = 1e-2 / problem.unc_num;
wc_samp_num = 1e6;
init_dual_points = init_dual_point';

[st1_deci, wc_sampler, primal_sol, dual_sol, DRO_UB, DRO_LB, error_sub, error_prob, output] ...
    = IterSolver.run(tolerance, init_dual_points, OT_tolerance);

RS = RandStream('mrg32k3a', 'Seed', 9999);
wc_sampler.setRandStream(RS);
wc_samples = wc_sampler.randSample(wc_samp_num);


fprintf('Total time = %.0fs\n', output.total_time);
fprintf('--- LSIP total time = %.0fs\n', output.LSIP_total_time);
fprintf('------- LSIP LP time = %.0fs\n', output.LSIP_LP_time);
fprintf('------- LSIP global time = %.0fs\n', output.LSIP_global_time);
fprintf('--- NewVar total time = %.0fs\n', output.NewVar_total_time);

save(CONFIG.SAVEPATH_OUTPUTS, ...
    'st1_deci', ...
    'wc_samples', ...
    'primal_sol', ...
    'dual_sol', ...
    'DRO_UB', ...
    'DRO_LB', ...
    'error_sub', ...
    'error_prob', ...
    'tolerance', ...
    'OT_tolerance', ...
    'output', ...
    '-v7.3');