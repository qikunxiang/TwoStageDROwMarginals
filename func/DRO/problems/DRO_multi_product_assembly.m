function DRO = DRO_multi_product_assembly(price_vec, assembly_mat, ...
    revenue_vec, salvage_vec, stage2_ub)
% Construct the DRO problem specification for a multi-product assembly
% problem
% Inputs: 
%       price_vec: vector containing the prices of parts
%       assembly_mat: matrix describing the parts needed for producing each
%       unit of product; each row represents a product and each column
%       represents a part
%       revenue_vec: vector containing the revenues from selling one unit
%       of each type of product
%       salvage_vec: vector containing the salvage prices of the remaining
%       parts
%       stage2_ub: specifies the upper bounds on the second-stage dual
%       variables
% Output:
%       DRO: struct containing the DRO specifications

K1 = length(price_vec);
N = length(revenue_vec);
K2ast = N + K1;

assert(length(salvage_vec) == K1);
assert(size(assembly_mat, 1) == N && size(assembly_mat, 2) == K1);
assert(length(stage2_ub) == K2ast);

DRO = struct;
DRO.stage1 = struct;
DRO.stage1.c1 = price_vec;
DRO.stage1.lb = zeros(K1, 1);
DRO.stage2 = struct;
DRO.stage2.V = [-speye(K1); sparse(N, K1)];
DRO.stage2.W = [sparse(K1, N); -speye(N)];
DRO.stage2.b = zeros(K2ast, 1);
DRO.stage2.Lin = [-sparse(assembly_mat), -speye(N)];
DRO.stage2.qin = -revenue_vec;
DRO.stage2.lb = [salvage_vec; zeros(N, 1)];
DRO.stage2.ub = stage2_ub;

end

