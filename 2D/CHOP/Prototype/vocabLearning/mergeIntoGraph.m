%> Name: mergeIntoGraph
%>
%> Description: Given the hierarchical graph and newly formed level, this
%> function embeds new level into the hierarchy by forming parent links
%> between current level and the previous level. If the graph structure has
%> position information, the positions of the children are averaged to find
%> out the position of their super-structure.
%>
%> @param graph Hierarchical graph structure.
%> @param level New level.
%> @param allLeafNodes Array containing leaf nodes.
%> @param levelItr Level iterator.
%> @param position Position calculations are processed if 1.
%>
%> @retval graph Newly formed hierarchical graph structure.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 02.12.2013
function [graph] = mergeIntoGraph(graph, level, allLeafNodes, levelItr, position)
    %% Go over children list of each instance in current level
    previousLevel = graph{levelItr-1};
    if position
        for newInstItr = 1:numel(level)
            leafNodes = level(newInstItr).leafNodes;
            % If necessary, process the position calculation/imageId copying work too.
            newPosition = [0,0];
            for childItr = 1:numel(leafNodes)
               newPosition = newPosition + allLeafNodes{leafNodes(childItr),2};
            end
            newPosition = round(newPosition / numel(leafNodes));
            level(newInstItr).position = newPosition;
        end
    end
    
    %% Reorder the nodes so that they are ordered by their image ids.
    if position
        imageIds = [level.imageId];
        [~, idx] = sort(imageIds);
        level = level(:,idx);
    end
    
    %% Assign the nodes in the previuos layer to their parents in this level.
    for newInstItr = 1:numel(level)
        children = level(newInstItr).children;
 %       leafNodes = [];
        for childItr = 1:numel(children)
           if ~ismember(newInstItr, previousLevel(children(childItr)).parents)
                previousLevel(children(childItr)).parents = ...
                    [previousLevel(children(childItr)).parents, newInstItr];
           end
        end 
    end
    
    %% Assign new levels and move on.
    graph(levelItr-1) = {previousLevel};
    graph(levelItr) = {level};
end