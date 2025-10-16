function CONFIG = SupplyNet_exp_config()
% Place to store global configurations of the two-stage distributionally robust supply chain network design problem with marginal
% constraints
% Output:
%   CONFIG: a struct containing configurations as fields

CONFIG = global_config();

% path to save files
CONFIG.SAVEPATH = [CONFIG.SAVEPATH_ROOT, 'DRO/SupplyNet_Exp/'];

% path to log files
CONFIG.LOGPATH = [CONFIG.LOGPATH_ROOT, 'DRO/SupplyNet_Exp/'];

% root folder for gurobi log files
CONFIG.LOGPATH_GUROBI = [CONFIG.LOGPATH, 'Gurobi/'];

% if the directory does not exist, create it first
if ~exist(CONFIG.SAVEPATH, 'dir')
    mkdir(CONFIG.SAVEPATH);
end

% if the directory does not exist, create it first
if ~exist(CONFIG.LOGPATH, 'dir')
    mkdir(CONFIG.LOGPATH);
end

% if the directory does not exist, create it first
if ~exist(CONFIG.LOGPATH_GUROBI, 'dir')
    mkdir(CONFIG.LOGPATH_GUROBI);
end

% file name of the inputs
CONFIG.FILENAME_INPUTS = 'inputs';
CONFIG.SAVEPATH_INPUTS = [CONFIG.SAVEPATH, CONFIG.FILENAME_INPUTS, '.mat'];

% file name of the outputs
CONFIG.FILENAME_OUTPUTS = 'outputs';
CONFIG.SAVEPATH_OUTPUTS = [CONFIG.SAVEPATH, CONFIG.FILENAME_OUTPUTS, '.mat'];

% file name of the main logs
CONFIG.LOGNAME_MAIN = 'main';
CONFIG.LOGPATH_MAIN = [CONFIG.LOGPATH, CONFIG.LOGNAME_MAIN, '.log'];

% file name of the iterative solver logs
CONFIG.LOGNAME_ITERSOLVER = 'iter_solver';
CONFIG.LOGPATH_ITERSOLVER = [CONFIG.LOGPATH, CONFIG.LOGNAME_ITERSOLVER, '.log'];

% file name of the linear semi-infinite programming logs
CONFIG.LOGNAME_LSIP_MAIN = 'LSIP_main';
CONFIG.LOGPATH_LSIP_MAIN = [CONFIG.LOGPATH, CONFIG.LOGNAME_LSIP_MAIN, '.log'];

% file name of the linear programming logs
CONFIG.LOGNAME_LSIP_LP = 'LSIP_LP';
CONFIG.LOGPATH_LSIP_LP = [CONFIG.LOGPATH_GUROBI, CONFIG.LOGNAME_LSIP_LP, '.log'];

% file name of the global optimization logs
CONFIG.LOGNAME_LSIP_GLOBAL = 'LSIP_global';
CONFIG.LOGPATH_LSIP_GLOBAL = [CONFIG.LOGPATH, CONFIG.LOGNAME_LSIP_GLOBAL, '.log'];

% file name of the new variable solver logs
CONFIG.LOGNAME_NEWVAR_MAIN = 'NewVar_main';
CONFIG.LOGPATH_NEWVAR_MAIN = [CONFIG.LOGPATH, CONFIG.LOGNAME_NEWVAR_MAIN, '.log'];

% file name of the new variable solver global logs
CONFIG.LOGNAME_NEWVAR_GLOBAL = 'NewVar_global';
CONFIG.LOGPATH_NEWVAR_GLOBAL = [CONFIG.LOGPATH_GUROBI, CONFIG.LOGNAME_NEWVAR_GLOBAL, '.log'];

end
