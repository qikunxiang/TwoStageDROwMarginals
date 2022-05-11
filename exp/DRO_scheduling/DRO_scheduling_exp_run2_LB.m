load('exp/DRO_scheduling/exp_inputs.mat');
load('exp/DRO_scheduling/exp_invcdf.mat');
load('exp/DRO_scheduling/exp_rst_UB.mat', ...
    'stage1_decision_cell', 'stage1_multipliers_cell', ...
    'coup_meas_cell');

% compute the inverse cdf of the marginal distributions
marginv = @(i, u)(interp1(mixtrnorm_invcdf{1}, mixtrnorm_invcdf{2}, u));

% Monte Carlo sample numbers and repetition numbers
samp_no = 1e7;
rep_no = 1000;
batch_size = 1e5;

LB_samps_cell = cell(step_no, 1);
LB_errbar = zeros(step_no, 2);
coup_meas_improved_cell = cell(step_no, 1);
stage1_multipliers_improved_cell = cell(step_no, 1);

for step_id = 1:step_no
    if isempty(coup_meas_cell{step_id})
        continue;
    end
    
    rand_stream = RandStream('combRecursive', 'Seed', 1000 + step_id);

    % collapse atoms that are identical in the x component
    [x_atoms, unique_atoms, prob_id] = unique( ...
        coup_meas_cell{step_id}.x, 'row');
    coup_improved = struct;
    coup_improved.x = x_atoms;
    coup_improved.probs = accumarray(prob_id, ...
        coup_meas_cell{step_id}.probs, [length(unique_atoms), 1]);

    [tighten_model, stage1_in_no, stage1_eq_no] ...
        = DRO_primal_tighten_gurobi(coup_improved.x, ...
        knots_cell{step_id}, coup_improved.probs, DRO);

    tighten_params = struct;
    tighten_params.OutputFlag = 1;

    tighten_output = gurobi(tighten_model, tighten_params);
    coup_improved.lambda = reshape(tighten_output.x( ...
        1:(K2ast * length(coup_improved.probs))), K2ast, ...
        length(coup_improved.probs))';
    stage1_multipliers_vec = tighten_output.x( ...
        K2ast * length(coup_improved.probs) + 1:end);
    stage1_multipliers_improved = struct;
    stage1_multipliers_improved.in ...
        = stage1_multipliers_vec(1:stage1_in_no);
    stage1_multipliers_improved.eq ...
        = stage1_multipliers_vec(stage1_in_no + 1:end);

    coup_meas_improved_cell{step_id} = coup_improved;
    stage1_multipliers_improved_cell{step_id} ...
        = stage1_multipliers_improved;

    fprintf('Setting %d:\n', step_id);
    LB_samps_cell{step_id} = DRO_primal_approx(DRO, ...
        stage1_multipliers_improved, coup_improved, ...
        marginv, samp_no, rep_no, batch_size, rand_stream, [], true);
    LB_errbar(step_id, :) = quantile( ...
        LB_samps_cell{step_id}, [0.025, 0.975]);
    
    save('exp/DRO_scheduling/exp_rst_LB.mat', ...
        'LB_samps_cell', 'LB_errbar', 'coup_meas_improved_cell', ...
        'stage1_multipliers_improved_cell', ...
        'samp_no', 'rep_no', 'batch_size', '-v7.3');
end