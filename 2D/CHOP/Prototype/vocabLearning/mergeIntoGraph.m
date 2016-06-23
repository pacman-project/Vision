%> Name: mergeIntoGraph
%>
%> Description: Given the hierarchical graph and newly formed level, this
%> function embeds new level into the hierarchy by forming parent links
%> between current level and the previous level. Only center nodes have
%> parent links: The reason being, we use parent info in the inference
%> procedure. Since center node is the only node that *has* to be present for
%> a part to be inferred, secondary nodes need not be linked to their
%> parents. Center node links are enough.
%>
%> @param graph Hierarchical graph structure.
%> @param level New level.
%> @param levelItr Level iterator.
%>
%> @retval graph Newly formed hierarchical graph structure.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 02.12.2013
function [graph] = mergeIntoGraph(graph, level, levelItr)
    %% Go over children list of each instance in current level
    previousLevel = graph{levelItr-1};
    
    %% Assign the nodes in the previuos layer to their parents in this level.
    emptyArr = {zeros(0, 'int32')};
    prevLevelParents = cell(numel(previousLevel),1);
    for itr = 1:numel(prevLevelParents)
        prevLevelParents(itr) = emptyArr;
    end
    previousLevelLabels = cat(1, previousLevel.label);
    for newInstItr = 1:numel(level)
        centerChild = level(newInstItr).children(1);
        % Find all low level nodes that correspond to this child.
        relevantPrevLevelNodes = find(previousLevelLabels == centerChild);
        for prevNodeItr = 1:numel(relevantPrevLevelNodes)
            curParents = prevLevelParents{relevantPrevLevelNodes(prevNodeItr)};
            if ~ismembc(int32(newInstItr), int32(curParents))
                curParents = cat(2, curParents, int32(newInstItr));
                prevLevelParents{relevantPrevLevelNodes(prevNodeItr)} = curParents;
            end
        end
    end
    [previousLevel.parents] = deal(prevLevelParents{:});
    
    %% Assign new levels and move on.
    graph(levelItr-1) = {previousLevel};
    graph(levelItr) = {level};
end