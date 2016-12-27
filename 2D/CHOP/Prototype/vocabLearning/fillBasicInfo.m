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
function graphLevel = fillBasicInfo(previousLevelImageIds, previousLevelLeafNodes, graphLevel, levelItr, options)
    numberOfNodes = numel(graphLevel);
    nodePositions = cat(1, graphLevel.precisePosition);
    
    % Learn stride.
    if strcmp(options.filterType, 'gabor')
         stride = options.gabor.stride;
    else
         stride = options.auto.stride;
    end
    poolDim = options.poolDim;
    
    % Calculate pool factor.
    poolFactor = nnz(~ismembc(2:levelItr, options.noPoolingLayers));
    clear options;
    
    childrenArr = {graphLevel.children};

    % Allocate space for fast data structures for current layer.
    newImageIds = cell(numberOfNodes,1);
    newPrecisePositions = newImageIds;
    newPositions = newImageIds;
    newLeafNodes = newImageIds;

    for newNodeItr = 1:numberOfNodes
        nodeChildren = childrenArr{newNodeItr};
        newImageIds{newNodeItr} = previousLevelImageIds(nodeChildren(1));
        nodeLeafNodes = cat(2, previousLevelLeafNodes{nodeChildren});
        nodeLeafNodes = fastsortedunique(sort(nodeLeafNodes));

        % Calculate both positions. For precise position, we obtain the
        % mean of the bounding box that is spanned by the leaf nodes.
  %      childrenPos = firstLevelPrecisePositions(nodeChildren, :);
  %      precisePosition = round((min(childrenPos,[], 1) + max(childrenPos, [], 1)) / 2);
  
        % We calculate new position based on leaf nodes.
%        precisePosition = round(mean(firstLevelPrecisePositions(nodeLeafNodes,:), 1));
  
        % Assign positions.
        precisePosition = nodePositions(newNodeItr, :);
        newPositions{newNodeItr} = int32(calculatePooledPositions(precisePosition, poolFactor, poolDim, stride));
        newLeafNodes{newNodeItr} = nodeLeafNodes;
    end

    % Finally, assign everything to the new structure.
    [graphLevel.imageId] = deal(newImageIds{:});
    [graphLevel.position] = deal(newPositions{:});
    [graphLevel.leafNodes] = deal(newLeafNodes{:});
    
    % Rearrange graph level so it is sorted by image id.
    arrayToSort = [[graphLevel.imageId]', [graphLevel.labelId]', cat(1, graphLevel.position)];
    [~, sortedIdx] = sortrows(arrayToSort);
    graphLevel = graphLevel(sortedIdx);
end