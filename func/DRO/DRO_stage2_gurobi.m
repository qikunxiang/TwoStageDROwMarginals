function model = DRO_stage2_gurobi(DRO)
% Process the constraints (including bounds) of the second-stage problem in
% a DRO problem and return a Gurobi model
% Inputs: 
%       DRO: the specification of the DRO problem
% Outputs: 
%       model: struct containing Gurobi-style constraints and bound
%       information

K2ast = size(DRO.stage2.V, 1);

% compute the truncation bounds for the conjugate functions
model = struct;

% generate the constraints for the second-stage decision
stage2_A = sparse(0, K2ast);
stage2_rhs = zeros(0, 1);
stage2_sense = [];

if isfield(DRO.stage2, 'Lin') && isfield(DRO.stage2, 'qin')
    % append the inequality constraints
    stage2_A = [stage2_A; sparse(DRO.stage2.Lin)];
    stage2_rhs = [stage2_rhs; DRO.stage2.qin];
    stage2_sense = [stage2_sense; repmat('<', length(DRO.stage2.qin), 1)];
end

if isfield(DRO.stage2, 'Leq') && isfield(DRO.stage2, 'qeq')
    % append the equality constraints
    stage2_A = [stage2_A; sparse(DRO.stage2.Leq)];
    stage2_rhs = [stage2_rhs; DRO.stage2.qeq];
    stage2_sense = [stage2_sense; repmat('=', length(DRO.stage2.qeq), 1)];
end

model.modelsense = 'max';
model.A = stage2_A;
model.rhs = stage2_rhs;
model.sense = stage2_sense;
model.lb = DRO.stage2.lb;
model.ub = DRO.stage2.ub;

end
