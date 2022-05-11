function DRO = DRO_task_scheduling(N, T)
% Construct the DRO problem specification for a task scheduling problem
% Inputs: 
%       N: the number of tasks
%       T: the length of the scheduling time window
% Output:
%       DRO: struct containing the DRO specifications

assert(T > 0);

DRO = struct;
DRO.stage1 = struct;
DRO.stage1.c1 = zeros(N, 1);
DRO.stage1.lb = zeros(N, 1);
DRO.stage1.Lin = sparse(ones(1, N));
DRO.stage1.qin = T;
DRO.stage2 = struct;
DRO.stage2.V = -speye(N);
DRO.stage2.W = speye(N);
DRO.stage2.b = zeros(N, 1);
DRO.stage2.Lin = sparse(repmat((1:N - 1)', 2, 1), [1:N - 1, 2:N]', ...
    [ones(N - 1, 1); -ones(N - 1, 1)], N - 1, N);
DRO.stage2.qin = ones(N - 1, 1);
DRO.stage2.lb = zeros(N, 1);
DRO.stage2.ub = (N:(-1):1)';

end

