classdef CGLP < ALGORITHM
% <2025> <multi> <real/integer/label/binary/permutation> <dynamic>
% Correlation-guided layered prediction
% FEinit --- 10000 --- Function evaluations for initialization
% taut   ---    10 --- Number of generations for static optimization

%------------------------------- Reference --------------------------------
% K. Yu, D. Zhang, J. Liang, K. Chen, C. Yue, K. Qiao, and L. Wang. A
% correlation-guided layered prediction approach for evolutionary dynamic
% multiobjective optimization. IEEE Transactions on Evolutionary
% Computation, 2025, 27(5): 1398-1412.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties
        POF_iter;
        POS_iter;
        PopX;
        Pareto;
    end
    methods
        function main(Algorithm, Problem)
            %% Parameter setting
            [FEinit,taut,restartFE] = Algorithm.ParameterSet(10000,10,0);
            % restartFE > 0：静态问题适配模式——每 restartFE 次评价强制唤醒一次
            %   （静态环境下 Changed 永不触发，周期性唤醒让预测机制参与优化，
            %    且 RMMEDA 的预算阈值改为递增目标，避免 FE 超过 FEinit 后空转）
            % restartFE = 0：原始动态模式（Changed 检测环境变化）
            nextWake  = FEinit;  % 静态模式：初始化完成后立即首次唤醒，之后按 restartFE 递增
            MaxT      = Problem.maxFE/taut;
            hisPareto = cell(MaxT,1);
            hisPopX   = {};
            hisPop    = {};
            tipe      = 0;
            % Reset the number of saved populations (only for dynamic optimization)
            Algorithm.save = sign(Algorithm.save)*inf;

            %% Generate random population
            T = 1;
            Population = Problem.Initialization();
            Algorithm.POF_iter = cell(1, Problem.maxFE/Problem.N);% 预分配内存
            Algorithm.POS_iter = cell(1, Problem.maxFE/Problem.N);
            AllPop = []; 
            for i = 1 : 2
                Population  = RMMEDA(Algorithm,Problem,Problem.N,FEinit);
                hisPopX{T}  = Algorithm.PopX';
                AllPop      = [AllPop,Population];
                hisPop{T}.F = Algorithm.Pareto.F';
                hisPop{T}.X = Algorithm.Pareto.X';
                [hisPareto,hisPop] = Get_C(hisPop,hisPareto,T);
                T = T + 1;
            end

            %% Optimization
            while Algorithm.NotTerminated(Population)
                if restartFE > 0
                    doWake = Problem.FE >= nextWake;      % 静态适配：周期性唤醒
                else
                    doWake = Changed(Problem,Population); % 原始动态模式：变化检测
                end
                if doWake
                    [Population1,pop_LCM,pop_DCM,tipe] = CGLP_pre(Problem,hisPop,T,hisPareto,Problem.N,tipe);
                    if restartFE > 0
                        targetFE = min(Problem.FE+restartFE,Problem.maxFE); % 递增预算目标
                    else
                        targetFE = FEinit;                                   % 原始逻辑
                    end
                    Population  = RMMEDA(Algorithm,Problem,Problem.N,targetFE,Population1);
                    tipe        = selfadjust(Algorithm.PopX',pop_LCM,pop_DCM,tipe);
                    hisPopX{T}  = Algorithm.PopX';
                    AllPop      = [AllPop,Population];
                    hisPop{T}.F = Algorithm.Pareto.F';
                    hisPop{T}.X = Algorithm.Pareto.X';
                    [hisPareto,hisPop] = Get_C(hisPop,hisPareto,T);
                    T = T+1;
                end
                if Problem.FE >= Problem.maxFE 
                    Population = [AllPop,Population]; 
                    [~,rank]=sort(Population.adds(zeros(length(Population),1))); 
                    Population = Population(rank); 
                end
             end                  
        end
    end  
end

function Population = RMMEDA(Algorithm, Problem,popSize,FEinit,init_pop)
    Lower = repmat(Problem.lower',1,popSize);
    Upper = repmat(Problem.upper',1,popSize);

    if nargin==4
        Population = Problem.Initialization();
    elseif nargin > 4
        if size(init_pop,2) < popSize
            newpop   = addNoise(init_pop, popSize, size(init_pop,2));
            init_pop = [init_pop' newpop'];
        else
            init_pop = init_pop';
        end
        init_pop   = max(min(init_pop,Upper),Lower);
        Population = Problem.Evaluation(init_pop');
    end

    iter = 1;
    while Problem.FE<FEinit
        Offspring  = Operator(Problem, Population);
        Population = EnvironmentalSelection([Population, Offspring], Problem.N);
        Algorithm.POF_iter{iter} = Population.objs;
        Algorithm.POS_iter{iter} = Population.decs;
        iter = iter + 1;
    end
    Algorithm.PopX     = Population.decs;
    Algorithm.Pareto.F = Population.objs;
    Algorithm.Pareto.X = Population.decs;
end