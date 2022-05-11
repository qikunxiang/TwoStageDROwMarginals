function [model, stage1_in_no, stage1_eq_no] ...
    = DRO_primal_tighten_gurobi(x_atoms, knots, probs, DRO)
% Formulate the problem of tightening the lower bound on a DRO problem
% given a discrete measure into an LP problem in Gurobi
% Inputs: 
%       x_atoms: the atoms of the discrete measure in the form of knot
%       indices (x part)
%       knots: cell array containing the knots in each dimension
%       probs: the corresponding probabilities
%       DRO: the specification of the DRO problem
% Outputs: 
%       model: the formulated problem in Gurobi
%       stage1_in_no: number of inequality constraints in the first-stage
%       problem
%       stage2_eq_no: number of equality constraints in the first-stage
%       problem

K1 = size(DRO.stage2.V, 2);

atom_no = size(x_atoms, 1);

assert(length(probs) == atom_no, 'the discrete measure is mis-specified');


knot_vec = vertcat(knots{:});
knot_cumnum = cumsum(cellfun(@length, knots));

% convert the x atoms from knot indices to actual locations
x_atoms = knot_vec(x_atoms + [0, knot_cumnum(1:end - 1)']);

% the objective vector
objvec = (DRO.stage2.W * x_atoms' + DRO.stage2.b) .* probs';

if isfield(DRO.stage1, 'Lin') && isfield(DRO.stage1, 'qin')
    objvec_in = DRO.stage1.qin;
    stage1_Ain = DRO.stage1.Lin;
    stage1_in_no = length(DRO.stage1.qin);
else
    objvec_in = zeros(0, 1);
    stage1_Ain = sparse(0, K1);
    stage1_in_no = 0;
end

if isfield(DRO.stage1, 'lb')
    % turn lower bounds into inequality constraints
    lb_constrained = ~isinf(DRO.stage1.lb);
    objvec_in = [objvec_in; -DRO.stage1.lb(lb_constrained)];
    stage1_Ain = [stage1_Ain; sparse(1:sum(lb_constrained), ...
        find(lb_constrained), -1, sum(lb_constrained), K1)];
    stage1_in_no = stage1_in_no + sum(lb_constrained);
end

if isfield(DRO.stage1, 'ub')
    % turn upper bounds into inequality constraints
    ub_constrained = ~isinf(DRO.stage1.ub);
    objvec_in = [objvec_in; DRO.stage1.ub(ub_constrained)];
    stage1_Ain = [stage1_Ain; sparse(1:sum(ub_constrained), ...
        find(ub_constrained), 1, sum(ub_constrained), K1)];
    stage1_in_no = stage1_in_no + sum(ub_constrained);
end

if isfield(DRO.stage1, 'Leq') && isfield(DRO.stage1, 'qeq')
    objvec_eq = DRO.stage1.qeq;
    stage1_Aeq = DRO.stage1.Leq;
    stage1_eq_no = length(DRO.stage1.qeq);
else
    objvec_eq = zeros(0, 1);
    stage1_Aeq = sparse(0, K1);
    stage1_eq_no = 0;
end

objvec = [objvec(:); objvec_in; objvec_eq];

% the constraints consist of two parts 

% the first part comes from the feasible set of the second-stage problem
stage2_model = DRO_stage2_gurobi(DRO);

A1_cell = cell(atom_no, 1);
A1_cell(1:end) = {stage2_model.A};

A1 = blkdiag(A1_cell{:});
A1 = [A1, sparse(size(A1, 1), stage1_in_no + stage1_eq_no)];
rhs1 = repmat(stage2_model.rhs, atom_no, 1);
sense1 = repmat(stage2_model.sense, atom_no, 1);

% the second part comes from the last set of constraints in the dual
% problem
A2_cell = cell(atom_no, 1);

for atomid = 1:atom_no
    A2_cell{atomid} = DRO.stage2.V' * probs(atomid);
end

A2 = [-horzcat(A2_cell{:}), stage1_Ain', stage1_Aeq'];
rhs2 = DRO.stage1.c1;
sense2 = repmat('=', K1, 1);

% lower and upper bounds
lb = [repmat(stage2_model.lb, atom_no, 1); ...
    -inf(stage1_in_no + stage1_eq_no, 1)];
ub = [repmat(stage2_model.ub, atom_no, 1); ...
    zeros(stage1_in_no, 1); inf(stage1_eq_no, 1)];

% assemble the model
model = struct;
model.modelsense = 'max';
model.obj = objvec;
model.A = [A1; A2];
model.rhs = [rhs1; rhs2];
model.sense = [sense1; sense2];
model.lb = lb;
model.ub = ub;

end