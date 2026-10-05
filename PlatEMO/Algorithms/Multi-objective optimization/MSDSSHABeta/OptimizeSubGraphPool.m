function SubGraphPool=OptimizeSubGraphPool(SubGraphPool,TempSubGraph,NumOfNonEmptySubgraphs)
%% Add Rule:
%       Rule1: Add subgraph if SubGraphPool has vacancies
%       Rule2: Randomly select a subgraph with fewer valid
%       bits(0,1) than current graph to replace
ValidBitCnt=sum(or(TempSubGraph(:)==0,TempSubGraph(:)==1));
if SubGraphPool.Limit- NumOfNonEmptySubgraphs>0
    AddIndex=NumOfNonEmptySubgraphs+1;
    SubGraphPool.Graphs{AddIndex,1}=TempSubGraph;
    SubGraphPool.ValidBitCnts(AddIndex,1)=ValidBitCnt;
    SubGraphPool.Size = SubGraphPool.Size+1;
else
%     AlternativeIndex=SubGraphPool.GetAlternativeIndex(ValidBitCnt);
    % index = find(SubGraphPool.ValidBitCnts(:)<ValidBitCnt);
    % if numel(index) ~=0
    %     AlternativeIndex=index(randperm(numel(index),1));
    % else
    %     AlternativeIndex=-1;   % 
    % end
    % if AlternativeIndex~=-1
    %     SubGraphPool.Graphs{AlternativeIndex,1}=TempSubGraph;
    %     SubGraphPool.ValidBitCnts(AlternativeIndex,1)=ValidBitCnt;        
    % end
    [~,index]=min(SubGraphPool.ValidBitCnts);
    sim_score=computeSimilarity(TempSubGraph,SubGraphPool.Graphs{index,1});
    if sim_score<0.5
        SubGraphPool.Graphs{index,1}=TempSubGraph;
        SubGraphPool.ValidBitCnts(index,1)=ValidBitCnt;
    end
end
end
function sim_score=computeSimilarity(H1,H2)
    num_matchs=sum(H1(:) == H2(:));
    total_elements = numel(H1);
    sim_score=num_matchs/total_elements;
end