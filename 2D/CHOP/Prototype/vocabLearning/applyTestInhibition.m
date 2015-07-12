%> Name: applyTestInhibition
%>
%> Description: Applies local inhibition of nodes to eliminate those who do
%> not introduce novelty, into the TEST graph. In this context, novelty is
%> measured by the ratio of leaf nodes not existing in the main graph over
%> all leaf nodes of the tested node. In case of an overlap, the node whose
%> label id is better wins. graphLevel is ASSUMED to be ordered by
%> imageIds, labelIds.
%>
%> @param graphLevel The current graph level, ASSUMED ordered by imageIds, 
%> then labelIds in an ascending manner.
%> @param options Program options.
%> @param levelItr Current level number.
%>
%> @retval graphLevel Modified graph level.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 11.12.2013
function [graphLevel] = applyTestInhibition(graphLevel, options, levelItr)
    % Calculate edge radius.
    scale = (1/options.scaling)^(levelItr-1);
    neighborhood = fix(options.edgeRadius * scale);
    noveltyThr = 1 - options.noveltyThr;
    downsampleRatio = floor((options.edgeQuantize-1)/2) / neighborhood;
    
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
    parfor imageId = 1:numberOfImages
        sortIdx = imageSortIdx{imageId};
        orgImageGraphLevel = imageGraphLevels{imageId};
        imageGraphLevel = orgImageGraphLevel(sortIdx);
        if isempty(imageGraphLevel)
            continue;
        end
        imageNodeCoords = imageAllNodeCoords{imageId};
        imageNodeLabels = cat(1, imageGraphLevel.labelId);
        numberOfNodesInImage = numel(imageGraphLevel);
        imagePreservedNodes = ones(numberOfNodesInImage,1)>0;
        imageLeafNodes = {imageGraphLevel.leafNodes}';
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
          distances = sqrt(sum(edgeCoords.^2, 2));
          adjacentNodes = imagePreservedNodes & distances <= neighborhood; 
          adjacentNodes(1:nodeItr) = 0;
          
          % Here, we do an elimination of nodes with the same id, that
          % happen to be very close to each other. It's mimicking the
          % pooling process in deep learning methods.
          normalizedEdgeCoords = fix(fix(downsampleRatio * double(edgeCoords(adjacentNodes, :))));
          eliminatedAdjacentNodes = normalizedEdgeCoords(:,1) == 0 & ...
                                        normalizedEdgeCoords(:,2) == 0 & ...
                                        imageNodeLabels(adjacentNodes) == imageNodeLabels(nodeItr);
          % If we need to eliminate some of the nodes based on spatial
          % adjacency, not because of overlapping leaf node support, we do
          % it here.
          if nnz(eliminatedAdjacentNodes) > 0
              adjacentNodeIdx = find(adjacentNodes);
              imagePreservedNodes(adjacentNodeIdx(eliminatedAdjacentNodes)) = 0;
              adjacentNodes(adjacentNodeIdx(eliminatedAdjacentNodes)) = 0;
          end
                                        
          % Go over each adjacent node, and apply inhibition if their leaf nodes are 
          % shared too much, under current novelty threshold.
          selfLeafNodes = imageLeafNodes{nodeItr};
          imagePreservedNodes(adjacentNodes) = cellfun(@(x,y) sum(ismembc(x, selfLeafNodes)) <= y, ...
              imageLeafNodes(adjacentNodes), maxSharedLeafNodes(adjacentNodes));
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
    clear imageGraphLevels imageAllNodeCoords preservedNodes;
end
