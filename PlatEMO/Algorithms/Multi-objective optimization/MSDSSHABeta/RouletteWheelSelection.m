function SubGraph = RouletteWheelSelection(SubGraphPool)
%ROULETTEWHEELSELECTION 在有效子图中进行轮盘赌选择
%   只在 ValidBitCnt > 0 的子图（即有固定位 0/1 的合法子图）中按固定位
%   数量占比做轮盘赌。先记录有效图的原始下标，抽中后再映射回 Graphs，
%   避免删除 0 概率项后压缩下标与池下标错位。

    idx = find(SubGraphPool.ValidBitCnts(1:SubGraphPool.Size) > 0);
    if isempty(idx)
        % 池中尚无合法子图时，退化为随机取一个已存在的图
        idx = (1:SubGraphPool.Size)';
    end
    p = SubGraphPool.ValidBitCnts(idx);
    if sum(p) == 0
        p = ones(numel(idx),1);
    end
    p = p ./ sum(p);
    SubGraph = SubGraphPool.Graphs{idx(find(rand <= cumsum(p), 1))};
end
