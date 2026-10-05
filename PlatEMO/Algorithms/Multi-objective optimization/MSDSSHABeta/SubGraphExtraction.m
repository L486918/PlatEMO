function [ZerosMap,OnesMap,SubGraphPool] = ...
                     SubGraphExtraction(SubGraphPool,~,ZerosMap,OnesMap,Problem,Population)
VarNum = sqrt(Problem.D);
TempSubGraph = zeros(VarNum);
IsExist= 0;

if Problem.FE <= Problem.maxFE/3
    alpha=0.6;
elseif Problem.FE > 2*Problem.maxFE/3
    alpha=0.7;
else
    alpha=0.9;
end

for i = 1 :VarNum
    for j = 1 : VarNum
        NumOfOnes=0;
        NumOfZeros=0;
        for  each = 1 : Problem.N
            if Population(each).dec(VarNum*(i-1)+j)==0
                NumOfZeros = NumOfZeros + 1;
            else
                NumOfOnes  = NumOfOnes + 1;
            end
        end
        ZerosMap(i,j) = NumOfZeros + ZerosMap(i,j);
        OnesMap(i,j) = NumOfOnes +  OnesMap(i,j);
        ProportionOfZeros = NumOfZeros/Problem.N;
        ProportionOfOnes = NumOfOnes/Problem.N;
        
        if ProportionOfOnes > alpha
            TempSubGraph(i,j)=1;
        elseif ProportionOfZeros > alpha
            TempSubGraph(i,j) = 0;
        else
            TempSubGraph(i,j) = 2;
        end
    end
end

% NumOfNonEmptySubgraphs = NumOfSubGraphs - ...
%         length(find(cell2mat(cellfun(@(x)(length(x)),SubGraphPool.Graphs,'Un',false))==0));
NumOfNonEmptySubgraphs=SubGraphPool.Size;
for each = 1:NumOfNonEmptySubgraphs
    if  isequal(TempSubGraph,SubGraphPool.Graphs{each,1})
        IsExist=1;
        break;
    end
end
if ~IsExist
    SubGraphPool=OptimizeSubGraphPool(SubGraphPool,TempSubGraph,NumOfNonEmptySubgraphs);
end


