classdef MSDSSHA < ALGORITHM
% <binary> 
% Main-correlation Statistics-based Decising Subgraph Heuristic Algorithm
% alpha --- 0.6 --- Parameter for extracting subgraph  

%------------------------------- Reference --------------------------------
% W. Li, X. Yao, T. Zhang, R. Wang, and L. Wang, Hierarchy ranking method
% for multimodal multi-objective optimization with local Pareto fronts,
% IEEE Transactions on Evolutionary Computation, 2022.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
methods
        function main(Algorithm, Problem)
            alpha = Algorithm.ParameterSet(0.6);
%             RootResultPath="C:\\Users\\Administrator\\Desktop\\result\\MSDSSHA";
%             RootAPath="C:\\code\\matlabCode\\MSDSS-HA\\ou\\";
%             subPath=readvars(fullfile(RootAPath,Problem.parameter.name,"a.txt"));
%             saveResultPath=fullfile(RootResultPath,num2str(subPath));
%             if exist(saveResultPath,'dir')==0
%                 mkdir(saveResultPath)
%             end
%             subPath=subPath+1;
%             modifyA(fullfile(RootAPath,Problem.parameter.name,"a.txt"),subPath,2);
           %% Generate random population
            Population          = Problem.Initialization();
            Offspring           = [];
            % Get the number of latent variables
            VarNum              = sqrt(Problem.D);
            % Set subgraph pool size
            SubGraphPool        = struct();
            SubGraphPool.Graphs = cell(100, 1);
            SubGraphPool.ValidBitCnts=zeros(100,1);
            SubGraphPool.Size=0;
            SubGraphPool.Limit=100;
            % Initialize maps for 0/1
            ZerosMap            = zeros(VarNum);
            OnesMap             = zeros(VarNum);

            %% Optimization
            % inFeasible 传 []：走基类缺省的 PEN 阈值统计；SubGraphPool 逐代存档
            while Algorithm.NotTerminated(Population,Offspring,[],SubGraphPool)
                [ZerosMap,OnesMap,SubGraphPool] = ...
                     SubGraphExtraction(SubGraphPool,alpha,ZerosMap,OnesMap,Problem,Population);
                 Offspring = ...
                        EnvironmentalSelection(alpha,OnesMap,ZerosMap,Population,Problem,SubGraphPool);
                Population=[Population,Offspring];
                [FrontNo,MaxFNo] = NDSort(Population.objs,Problem.N); % 2n个个体
                Next = FrontNo < MaxFNo;
                Last=find(FrontNo(:)==MaxFNo);
                % 决策空间拥挤度截断替代随机选择：保留结构上更稀疏的个体，提升多样性
                CrowdDis = CrowdingDistance(Population.decs,FrontNo);
                [~,Rank] = sort(CrowdDis(Last),'descend');
                Next(Last(Rank(1:Problem.N-sum(Next)))) = true;
                Population = Population(Next);
            end
        end
    end
end



