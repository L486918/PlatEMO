function Offspring = ...
        EnvironmentalSelection(alpha,OnesMap,ZerosMap,Population,Problem,SubGraphPool)
    if mod(Problem.FE/Problem.N,2) == 0
        Offspring = OperatorGA(Problem,Population);
    else
        beta= (Problem.FE*log10(exp(1))) / alpha;
        Offspring=zeros(Problem.N,Problem.D);
        % 阈值掩码每代只算一次：未决位上与统计趋势一致则强制填充
        Mask1 = OnesMap>=beta & ZerosMap<beta;   % 强制填1
        Mask0 = ZerosMap>=beta & OnesMap<beta;   % 强制填0
        for t = 1 : Problem.N
            subGraph = RouletteWheelSelection(SubGraphPool);
            % 未决位(==2)批量填充：先按掩码固化，其余随机（向量化）
            undecided = subGraph==2;
            subGraph(undecided & Mask1) = 1;
            subGraph(undecided & Mask0) = 0;
            rest = undecided & ~Mask1 & ~Mask0;
            subGraph(rest) = randi([0 1], nnz(rest), 1);
            Offspring(t,:)=reshape(subGraph',1,[]);
        end
        Offspring = Problem.Evaluation(Offspring);
    end
    
end

