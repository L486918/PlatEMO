function Offspring = ...
        EnvironmentalSelection(alpha,OnesMap,ZerosMap,Population,Problem,SubGraphPool)
    if mod(Problem.FE/Problem.N,2) == 0
        Offspring = OperatorGA(Problem,Population);
    else
        beta= (Problem.FE*log10(exp(1))) / alpha;
        Offspring=zeros(Problem.N,Problem.D);
        maxIter = Problem.maxFE / Problem.N; % 假设Problem定义了最大迭代次数
        currentIter = Problem.FE / Problem.N;
        if and(currentIter >= 0.3 * maxIter, SubGraphPool.Size>=3)

            selectionStrategy = 'tournament'; % 后期使用锦标赛选择
        else
            selectionStrategy = 'fitness_sharing'; % 早期使用适应度共享
        end
        % 阈值掩码每代只算一次：未决位上与统计趋势一致则强制填充
        Mask1 = OnesMap>=beta & ZerosMap<beta;   % 强制填1
        Mask0 = ZerosMap>=beta & OnesMap<beta;   % 强制填0
        for t = 1 : Problem.N
            % 动态选择策略
            if strcmp(selectionStrategy, 'fitness_sharing')
                subGraph = FitnessSharingSelection(SubGraphPool, Population, OnesMap, ZerosMap, beta);
            else
                    subGraph = TournamentSelection(SubGraphPool, 3); % 锦标赛大小为3
            end

            %subGraph = RouletteWheelSelection(SubGraphPool);
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
function selected = FitnessSharingSelection(pool, population, OnesMap, ZerosMap, beta)
    % 计算适应度共享后的概率分布
    sharedFitness = calculateSharedFitness(pool, population, OnesMap, ZerosMap, beta);
    
    % 根据共享后的适应度选择子图
    prob = sharedFitness / sum(sharedFitness);
    selectedIdx = randsample(1:pool.Size, 1, true, prob);
    selected = pool.Graphs{selectedIdx,1};
end
function sharedFitness = calculateSharedFitness(pool, ~, OnesMap, ZerosMap, beta)
    % 计算每个子图的原始适应度
    rawFitness=zeros(1,pool.Size);
    for i = 1:pool.Size
        rawFitness(i) = evaluateSubGraphFitness(pool.Graphs{i,1}, OnesMap, ZerosMap, beta);
    end
    % 计算相似性并调整适应度
    % 距离用归一化汉明距离（不同位占比，取值[0,1]），sigma=0.2 才有实际
    % 共享作用；原来的欧氏距离在 0/1/2 矩阵上恒 >= 1，惩罚永不触发。
    sigma = 0.2; % 共享半径（归一化汉明距离口径，建议范围 0.1~0.3）
    sharedFitness = rawFitness;
    for i = 1:pool.Size
        for j = 1:pool.Size
            if i ~= j
                distance = mean(pool.Graphs{i}(:) ~= pool.Graphs{j}(:)); % 归一化汉明距离
                if distance < sigma
                    sharedFitness(i) = sharedFitness(i) * (1 - distance/sigma);
                end
            end
        end
    end
end


function fitness = evaluateSubGraphFitness(subGraph, OnesMap, ZerosMap, beta)
    % 评估子图质量：统计符合阈值条件的变量数量
    validOnes = sum(subGraph(:) == 1 & OnesMap(:) > beta);
    validZeros = sum(subGraph(:) == 0 & ZerosMap(:) > beta);
    fitness = validOnes + validZeros;
end

%% 辅助函数：锦标赛选择
function selected = TournamentSelection(pool, tournamentSize)
    candidates = randperm(pool.Size, tournamentSize);
    fitnessValues = cellfun(@(x) sum(x(:)), pool.Graphs(candidates));
    [~, idx] = max(fitnessValues);  %% 这里选择的就是最浅的子图
    selected = pool.Graphs{candidates(idx)};
end




