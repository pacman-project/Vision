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

    % Fill in necessary internal structures.
    imageIds = [graphLevel.imageId];
    numberOfImages = max(imageIds);
    
    % Get node coordinates.
    nodeCoords = cat(1, graphLevel.position);
    
    %% Prepare data structures for parallel processing.
    imageGraphLevels = cell(numberOfImages,1);
    imageAllNodeCoords = cell(numberOfImages,1);
    for imageId = 1:numberOfImages
        imageNodeIdx = imageIds == imageId;
        imageGraphLevels(imageId) = {graphLevel(1,imageNodeIdx)};
        imageAllNodeCoords(imageId) = {nodeCoords(imageNodeIdx,:)}; 
    end
    
    %% Go over each node and check neighboring nodes for novelty introduced. Eliminate weak ones.
    preservedNodes = cell(numberOfImages,1);
    parfor imageId = 1:numberOfImages
        imageGraphLevel = imageGraphLevels{imageId};
        imageNodeCoords = imageAllNodeCoords{imageId};
        numberOfNodesInImage = numel(imageGraphLevel);
        imagePreservedNodes = ones(numberOfNodesInImage,1)>0;
        imageLeafNodes = {imageGraphLevel.leafNodes}';
        maxSharedLeafNodes = cellfun(@(x) numel(x) * noveltyThr , imageLeafNodes, 'UniformOutput', false);
        
        for nodeItr = 1:(numberOfNodesInImage-1)
          %% If nobody has erased this node before, it has a right to be in the final graph.
          if imagePreservedNodes(nodeItr) == 0
              continue;
          end
          
          %% Get each neighboring node.
          thisNodeCoords = imageNodeCoords(nodeItr,:);
          centerArr = repmat(thisNodeCoords, numberOfNodesInImage, 1);
          distances = sqrt(sum((centerArr - imageNodeCoords).^2, 2));
          adjacentNodes = imagePreservedNodes & distances <= neighborhood; 
          adjacentNodes(1:nodeItr) = 0;
          selfLeafNodes = imageLeafNodes{nodeItr};
          
          %% Go over each adjacent node, and apply inhibition if their leaf nodes are too common under current novelty threshold.
          imagePreservedNodes(adjacentNodes) = cellfun(@(x,y) sum(ismembc(x, selfLeafNodes)) <= y, ...
              imageLeafNodes(adjacentNodes), maxSharedLeafNodes(adjacentNodes));
        end
        preservedNodes(imageId) = {imagePreservedNodes'};
    end
    preservedNodes = [preservedNodes{:}]>0;
    graphLevel = graphLevel(:,preservedNodes);
end
