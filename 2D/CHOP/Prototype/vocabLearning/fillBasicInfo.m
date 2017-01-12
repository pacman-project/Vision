%> Name: fillBasicInfo
%>
%> Description: Given current level and previous level (both object graphs), 
%> fill in the current level's necessary basic info, including image id,
%> position, and leafNodes. In the end, the output is sorted based on imageId
%> and labelId.
%>
%> @param previousLevel 
%> @param graphLevel 
%> @param levelItr 
%> @param numberOfThreads
%>  
%> @retval graphLevel 
%>
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 04.02.2014
function graphLevel = fillBasicInfo(previousLevelImageIds, previousLevelLeafNodes, previousLevelPrecisePositions, graphLevel)
    numberOfNodes = numel(graphLevel);
    clear options;
    
    childrenArr = {graphLevel.children};

    % Allocate space for fast data structures for current layer.
    newImageIds = cell(numberOfNodes,1);
    newPositions = newImageIds;
    newLeafNodes = newImageIds;

    for newNodeItr = 1:numberOfNodes
        nodeChildren = childrenArr{newNodeItr};
        newImageIds{newNodeItr} = previousLevelImageIds(nodeChildren(1));
        nodeLeafNodes = cat(2, previousLevelLeafNodes{nodeChildren});
        nodeLeafNodes = fastsortedunique(sort(nodeLeafNodes));

        % Assign positions.
        childPositions = previousLevelPrecisePositions(nodeChildren, :);
        newPositions{newNodeItr} = round((min(childPositions) + max(childPositions))/2);
        newLeafNodes{newNodeItr} = nodeLeafNodes;
    end

    % Finally, assign everything to the new structure.
    [graphLevel.imageId] = deal(newImageIds{:});
    [graphLevel.precisePosition] = deal(newPositions{:});
    [graphLevel.leafNodes] = deal(newLeafNodes{:});
    
    % Rearrange graph level so it is sorted by image id.
    arrayToSort = [[graphLevel.imageId]', [graphLevel.labelId]', cat(1, graphLevel.position)];
    [~, sortedIdx] = sortrows(arrayToSort);
    graphLevel = graphLevel(sortedIdx);
end