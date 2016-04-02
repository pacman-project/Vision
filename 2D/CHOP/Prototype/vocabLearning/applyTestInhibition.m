%> Name: applyTestInhibition
%>
%> Description: Applies local inhibition of nodes to eliminate those who do
%> not introduce novelty, into the TEST graph. In this context, novelty is
%> measured by the ratio of leaf nodes not existing in the main graph over
%> all leaf nodes of the tested node. In case of an overlap, the node whose
%> label id is smaller wins. graphLevel is ASSUMED to be ordered by
%> imageIds, importance, labelIds. The importance of every graph is usually
%> linked to its activation. 
%>
%> @param graphLevel The current graph level, ASSUMED ordered by imageIds, 
%> then activations in a descending manner.
%> @param options Program options.
%> @param levelItr Current level id.
%>
%> @retval graphLevel Remaining graph level nodes.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 11.12.2013
function [graphLevel] = applyTestInhibition(graphLevel, options)
    % Calculate edge radius.
    noveltyThr = 1 - options.noveltyThr;
    fullOverlapThr = 0.99;
    halfMatrixSize = (options.receptiveFieldSize-1)/2 - 2;
    
    if numel(graphLevel) == 0
        return;
    end
    
    % Fill in necessary internal structures.
    imageIds = [graphLevel.imageId];
    numberOfImages = max(imageIds);
    
    if isempty(numberOfImages) || numberOfImages < 1
       return; 
    end
    
    % Get node coordinates.
    nodeCoords = cat(1, graphLevel.position);
    
    %% Prepare data structures for parallel processing.
    imageGraphLevels = cell(numberOfImages,1);
    imageSortIdx = cell(numberOfImages,1);
    imageAllNodeCoords = cell(numberOfImages,1);
    for imageId = 1:numberOfImages
        imageNodeIdx = imageIds == imageId;
        imageNodeCoords = nodeCoords(imageNodeIdx,:);
        if nnz(imageNodeIdx) == 0
            continue;
        end
        imageGraphLevel = graphLevel(1,imageNodeIdx);
        
        % Sort nodes based on their activations, and preserve the ordering.
        activationArr = [imageGraphLevel.activation];
        [~, sortIdx] = sort(activationArr, 'descend');
        imageGraphLevels(imageId) = {imageGraphLevel};
        imageSortIdx(imageId) = {sortIdx};
        
        % Save image coords.
        imageAllNodeCoords(imageId) = {imageNodeCoords(sortIdx,:)}; 
    end
    
    %% Go over each node and check neighboring nodes for novelty introduced. Eliminate weak ones.
    for imageId = 1:numberOfImages
        sortIdx = imageSortIdx{imageId};
        orgImageGraphLevel = imageGraphLevels{imageId};
        imageGraphLevel = orgImageGraphLevel(sortIdx);
        if isempty(imageGraphLevel)
            continue;
        end
        imageNodeCoords = imageAllNodeCoords{imageId};
        numberOfNodesInImage = numel(imageGraphLevel);
        imagePreservedNodes = ones(numberOfNodesInImage,1)>0;
        imageLeafNodes = {imageGraphLevel.leafNodes}';
        imageLeafNodeCounts = cellfun(@(x) numel(x), imageLeafNodes);
        maxSharedLeafNodes = cellfun(@(x) numel(x) * noveltyThr , imageLeafNodes, 'UniformOutput', false);
        
        for nodeItr = 1:(numberOfNodesInImage-1)
          %% If nobody has erased this node before, it has a right to be in the final graph.
          if imagePreservedNodes(nodeItr) == 0
              continue;
          end
          
          %% Get each neighboring node's relative coords.
          thisNodeCoords = imageNodeCoords(nodeItr,:);
          centerArr = repmat(thisNodeCoords, numberOfNodesInImage, 1);
          edgeCoords = imageNodeCoords - centerArr;
          
          % Find distances relative to the center node, and then obtain
          % leaf node supports.
          adjacentNodes = imagePreservedNodes & ... 
                edgeCoords(:,1) > -halfMatrixSize & ...
                edgeCoords(:,1) < halfMatrixSize & ...
                edgeCoords(:,2) > -halfMatrixSize & ...
                edgeCoords(:,2) < halfMatrixSize;
          adjacentNodes(1:nodeItr) = 0;
              
          % Go over each adjacent node, and apply inhibition if their leaf nodes are 
          % shared too much, under current novelty threshold.
          selfLeafNodes = imageLeafNodes{nodeItr};
          if noveltyThr >= fullOverlapThr
               % As an additional check, if the overlap threshold is really
               % small, we only eliminate fully overlapping nodes, that is, 
               % both overlapping with one another %100.
               imagePreservedNodes(adjacentNodes) = cellfun(@(x,y) sum(ismembc(x, selfLeafNodes)) <= y, ...
                   imageLeafNodes(adjacentNodes), maxSharedLeafNodes(adjacentNodes)) | ...
                   imageLeafNodeCounts(adjacentNodes) ~= numel(selfLeafNodes);
          else
              imagePreservedNodes(adjacentNodes) = cellfun(@(x,y) sum(ismembc(x, selfLeafNodes)) <= y, ...
                   imageLeafNodes(adjacentNodes), maxSharedLeafNodes(adjacentNodes));
          end
        end
        sortIdx = sort(sortIdx(imagePreservedNodes));
        imageGraphLevel = orgImageGraphLevel(sortIdx);
        imageGraphLevels{imageId} = imageGraphLevel;
    end
    graphLevel = [imageGraphLevels{:}];
    
    % Rearrange graph level so it is sorted by image id.
    arrayToSort = [[graphLevel.imageId]', [graphLevel.labelId]'];
    [~, sortedIdx] = sortrows(arrayToSort);
    graphLevel = graphLevel(sortedIdx);
    clearvars -except graphLevel
end
