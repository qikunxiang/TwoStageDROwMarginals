load('exp/DRO_scheduling/exp_inputs.mat');

LSIP_LB_list = zeros(step_no, 1);
LSIP_UB_list = zeros(step_no, 1);

dual_func_cell = cell(step_no, 1);
coup_meas_cell = cell(step_no, 1);
stage1_decision_cell = cell(step_no, 1);
stage1_multipliers_cell = cell(step_no, 1);
output_cell = cell(step_no, 1);
options_cell = cell(step_no, 1);

options = struct;
options.tolerance = LSIP_tolerance;
options.display = true;
options.log_file = 'scheduling_cutplane.log';

options.reduce_cuts = struct;
options.reduce_cuts.thres = 200;
options.reduce_cuts.freq = 200;

options.LP_params = struct;
options.LP_params.OutputFlag = 1;
options.LP_params.LogToConsole = 0;
options.LP_params.LogFile = 'scheduling_gurobi_LP.log';

options.global_params = struct;
options.global_params.OutputFlag = 1;
options.global_params.LogToConsole = 0;
options.global_params.LogFile = 'scheduling_gurobi_MILP.log';

options.rescue = struct;
options.rescue.save_file = 'scheduling_savefile.mat';
options.rescue.save_interval = 200;

reduce_cuts_thres = [200; 50; 50; 50; 20; 20; 10; 10];


for step_id = 1:step_no
    if step_id == 1
        coup_init = coup_init1;
    else
        knots_ind = cell(N, 1);

        for i = 1:N
            knots_ind{i} = (1:length(knots_cell{step_id}{i}))';
        end

        prev_atoms = coup_meas_cell{step_id - 1}.x;
        prev_probs = coup_meas_cell{step_id - 1}.probs;
        
        % collapse atoms that are identical in the x component
        [prev_atoms, unique_atoms, prob_id] = unique(prev_atoms, 'row');
        prev_probs = accumarray(prob_id, prev_probs, ...
            [length(unique_atoms), 1]);


        [atoms_new, probs_new] = reassembly_increment( ...
            knots_ind, expval_cell{step_id}, ...
            prev_atoms, prev_probs);
        coup_init = struct;
        [coup_init.x, unique_atoms, prob_id] = unique(atoms_new, 'row');
        coup_init.probs = accumarray(prob_id, probs_new, ...
            [length(unique_atoms), 1]);
    end


    tighten_model = DRO_primal_tighten_gurobi(coup_init.x, ...
        knots_cell{step_id}, coup_init.probs, DRO);

    tighten_params = struct;
    tighten_params.OutputFlag = 1;

    tighten_output = gurobi(tighten_model, tighten_params);
    coup_init.lambda = reshape(tighten_output.x( ...
        1:(K2ast * length(coup_init.probs))), K2ast, ...
        length(coup_init.probs))';

    options.reduce_cuts.thres = reduce_cuts_thres(step_id);

    [LSIP_UB_list(step_id), LSIP_LB_list(step_id), ...
        coup_meas_cell{step_id}, stage1_decision_cell{step_id}, ...
        stage1_multipliers_cell{step_id}, dual_func_cell{step_id}, ...
        output_cell{step_id}, options_cell{step_id}] ...
        = DRO_cutplane(DRO, knots_cell{step_id}, expval_cell{step_id}, ...
        coup_init, options);

    save('exp/DRO_scheduling/exp_rst_UB.mat', 'LSIP_LB_list', ...
        'LSIP_UB_list', 'dual_func_cell', 'coup_meas_cell', ...
        'stage1_decision_cell', 'stage1_multipliers_cell', ...
        'output_cell', 'options_cell', '-v7.3');
end