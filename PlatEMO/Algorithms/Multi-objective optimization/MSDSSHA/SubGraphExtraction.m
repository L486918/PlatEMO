function [ZerosMap,OnesMap,SubGraphPool] = ...
                     SubGraphExtraction(SubGraphPool,alpha,ZerosMap,OnesMap,Problem,Population)
% 论文 Eq.26-28 的忠实实现：
%   ZerosMap/OnesMap 为 Ξ0/Ξ1 的累计计数矩阵（跨代累加）；
%   子图提取使用累计比例 φk(r) = 累计计数 / 累计试验数(FE ≈ N·r)。

VarNum = sqrt(Problem.D);

% 当前代每位取 1 的次数（向量化，替代原三重循环）
% 注意 dec 的下标 VarNum*(i-1)+j 是行主序，reshape 后需转置对齐
OnesCnt  = reshape(sum(Population.decs,1), VarNum, VarNum).';
ZerosCnt = Problem.N - OnesCnt;

% 累计更新 Ξ1/Ξ0（Eq.26-27）
OnesMap  = OnesMap  + OnesCnt;
ZerosMap = ZerosMap + ZerosCnt;

% 累计比例 φ1(r)、φ0(r)（Eq.27：前 r 次迭代全部 N·r 个样本中的占比）
ProportionOfOnes  = OnesMap  / Problem.FE;
ProportionOfZeros = ZerosMap / Problem.FE;

% 指示函数生成子图（Eq.28）：1 优先，其次 0，否则未决位 2
TempSubGraph = 2*ones(VarNum);
TempSubGraph(ProportionOfZeros > alpha) = 0;
TempSubGraph(ProportionOfOnes  > alpha) = 1;

% 去重后加入子图池
IsExist = false;
for each = 1 : SubGraphPool.Size
    if isequal(TempSubGraph,SubGraphPool.Graphs{each,1})
        IsExist = true;
        break;
    end
end
if ~IsExist
    SubGraphPool = OptimizeSubGraphPool(SubGraphPool,TempSubGraph,SubGraphPool.Size);
end
end
