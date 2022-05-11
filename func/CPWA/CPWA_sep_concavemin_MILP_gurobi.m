function model = CPWA_sep_concavemin_MILP_gurobi(const, lin_vec, ...
    trans_mat, knots, knots_val, constraints)
% Formulate a separable concave piece-wise affine minimization problem into
% a mixed-integer linear programming problem in Gurobi
% Inputs: 
%       const: the constant intercept
%       lin_vec: the linear part of the objective function in a vector
%       trans_mat: the transformation matrix before applying the separable
%       concave function
%       knots: cell array containing the knots for the one-dimensional CPWA
%       functions in each dimension
%       knots_val: cell array containing the values of the one-dimensional
%       CPWA functions in each dimension
%       constraints: struct with the following fields
%           Lin: the inequality constraint matrix (<=) 
%           qin: the inequality constraint rhs
%           Leq: the equality constraint matrix
%           qeq: the equality constraint rhs
%           lb: the variable lower bounds
%           ub: the variable upper bounds
% Outputs: 
%       model: the formulated model in Gurobi

K2ast = length(lin_vec);
N = length(knots);
knotno_list = cellfun(@length, knots);
knotno = sum(knotno_list);

% check the inputs
assert(size(trans_mat, 1) == N && size(trans_mat, 2) == K2ast, ...
    'the transformation matrix mis-specified');
assert(all(cellfun(@length, knots_val) == knotno_list), ...
    'the values of one-dimensional functions mis-specified');

assert(length(constraints.lb) == K2ast);
assert(length(constraints.ub) == K2ast);

% decision variables summary:
%       name        type                length
%       ------------  Linear part ------------------
%       lambda      continuous          K2ast
%       ------------  Separable part  --------------
%       w           continuous          N
%       z_ij        continuous          knotno - N
%       iota_ij     binary              knotno - 2 * N

% compute the length of the decision variable
varlength = K2ast + N + knotno - N + knotno - 2 * N;

objcon = const;
objvec = zeros(varlength, 1);
objvec(1:K2ast) = lin_vec;

lb = -inf(varlength, 1);
lb(1:K2ast) = constraints.lb;

ub = inf(varlength, 1);
ub(1:K2ast) = constraints.ub;

vtype = [repmat('C', K2ast, 1); repmat('C', N, 1);
    repmat('C', knotno - N, 1); repmat('B', knotno - 2 * N, 1)];

% the linear part

% constraints summary
%   description                         number of constraints
%   lambda inequality constraints       length(constraints.qin)
%   lambda equality constraints         length(constraints.qin)
%   equalities relating lambda and w    N

lambda_A = sparse(0, K2ast);
lambda_rhs = zeros(0, 1);
lambda_sense = [];

if isfield(constraints, 'Lin') && isfield(constraints, 'qin')
    assert(size(constraints.Lin, 2) == K2ast);
    assert(size(constraints.Lin, 1) == length(constraints.qin));
    lambda_A = [lambda_A; sparse(constraints.Lin)];
    lambda_rhs = [lambda_rhs; constraints.qin];
    lambda_sense = [lambda_sense; repmat('<', length(constraints.qin), 1)];
end

if isfield(constraints, 'Leq') && isfield(constraints, 'qeq')
    assert(size(constraints.Leq, 2) == K2ast);
    assert(size(constraints.Leq, 1) == length(constraints.qeq));
    lambda_A = [lambda_A; sparse(constraints.Leq)];
    lambda_rhs = [lambda_rhs; constraints.qeq];
    lambda_sense = [lambda_sense; repmat('=', length(constraints.qeq), 1)];
end

lambda_A = [lambda_A, sparse(length(lambda_rhs), varlength - K2ast)];

link_A = [sparse(trans_mat), -speye(N), sparse(N, varlength - K2ast - N)];
link_rhs = zeros(N, 1);
link_sense = repmat('=', N, 1);

% the separable part

% constraints summary
%   description             number of constraints
%   relating w and z        N
%   relating z and iota     2 * (knotno - 2 * N)

sep1_rowno = N;
sep1_entryno = knotno;
sep1_r = zeros(sep1_entryno, 1);
sep1_c = zeros(sep1_entryno, 1);
sep1_v = zeros(sep1_entryno, 1);
sep1_rhs = zeros(sep1_rowno, 1);
sep1_sense = repmat('=', sep1_rowno, 1);

sep2_rowno = 2 * (knotno - 2 * N);
sep2_entryno = 2 * sep2_rowno;
sep2_r = zeros(sep2_entryno, 1);
sep2_c = zeros(sep2_entryno, 1);
sep2_v = zeros(sep2_entryno, 1);
sep2_rhs = zeros(sep2_rowno, 1);
sep2_sense = repmat('<', sep2_rowno, 1);

zcounter = K2ast + N;
iotacounter = K2ast + N + sum(knotno_list) - N;
sep1counter = 0;
sep2counter = 0;
sep2rowcounter = 0;

for i = 1:N
    lb(K2ast + i) = knots{i}(1);
    ub(K2ast + i) = knots{i}(end);
    
    % set objective
    objcon = objcon + knots_val{i}(1);
    objvec(zcounter + (1:knotno_list(i) - 1)) = diff(knots_val{i});
    
    % add bounds
    ub(zcounter + 1) = 1;
    lb(zcounter + knotno_list(i) - 1) = 0;
    
    % set the constraints relating x and z
    sep1_r(sep1counter + (1:knotno_list(i))) = i;
    sep1_c(sep1counter + 1) = K2ast + i;
    sep1_v(sep1counter + 1) = 1;
    sep1_c(sep1counter + (2:knotno_list(i))) = ...
        zcounter + (1:knotno_list(i) - 1);
    sep1_v(sep1counter + (2:knotno_list(i))) = -diff(knots{i});
    sep1_rhs(i) = knots{i}(1);
    
    % set the constraints relating z and iota
    sep2_r(sep2counter + (1:4 * (knotno_list(i) - 2))) = ...
        sep2rowcounter + repelem((1:2 * (knotno_list(i) - 2))', 2, 1);
    sep2_c(sep2counter + (1:4:4 * (knotno_list(i) - 2))) = ...
        zcounter + (1:knotno_list(i) - 2);
    sep2_v(sep2counter + (1:4:4 * (knotno_list(i) - 2))) = -1;
    sep2_c(sep2counter + (2:4:4 * (knotno_list(i) - 2))) = ...
        iotacounter + (1:knotno_list(i) - 2);
    sep2_v(sep2counter + (2:4:4 * (knotno_list(i) - 2))) = 1;
    sep2_c(sep2counter + (3:4:4 * (knotno_list(i) - 2))) = ...
        zcounter + (2:knotno_list(i) - 1);
    sep2_v(sep2counter + (3:4:4 * (knotno_list(i) - 2))) = 1;
    sep2_c(sep2counter + (4:4:4 * (knotno_list(i) - 2))) = ...
        iotacounter + (1:knotno_list(i) - 2);
    sep2_v(sep2counter + (4:4:4 * (knotno_list(i) - 2))) = -1;
    sep2_rhs = zeros(sep2_rowno, 1);
    
    
    % update the counters
    zcounter = zcounter + knotno_list(i) - 1;
    iotacounter = iotacounter + knotno_list(i) - 2;
    sep1counter = sep1counter + knotno_list(i);
    sep2counter = sep2counter + (knotno_list(i) - 2) * 4;
    sep2rowcounter = sep2rowcounter + (knotno_list(i) - 2) * 2;
end

sep1_A = sparse(sep1_r, sep1_c, sep1_v, sep1_rowno, varlength);
sep2_A = sparse(sep2_r, sep2_c, sep2_v, sep2_rowno, varlength);
sep_A = [sep1_A; sep2_A];
sep_rhs = [sep1_rhs; sep2_rhs];
sep_sense = [sep1_sense; sep2_sense];

% assemble the Gurobi model
model = struct;

model.modelsense = 'min';
model.objcon = objcon;
model.obj = objvec;
model.lb = lb;
model.ub = ub;
model.vtype = vtype;
model.A = [lambda_A; link_A; sep_A];
model.rhs = [lambda_rhs; link_rhs; sep_rhs];
model.sense = [lambda_sense; link_sense; sep_sense];

end
