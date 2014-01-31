%> Name: applyLocalInhibition
%>
%> Description: Applies local inhibition of nodes to eliminate those who do
%> not introduce novelty, into the graph. In this context, novelty is
%> measured by the ratio of leaf nodes not existing in the main graph over
%> all leaf nodes of the tested node. In case of an overlap, the node whose
%> label id is better wins. graphLevel is ASSUMED to be ordered by
%> labelIds.
%>
%> @param graphLevel The current graph level, ASSUMED ordered by labelIds.
%> @param options Program options.
%> @param levelItr Current level number.
%>
%> @retval graphLevel Modified graph level.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 11.12.2013
function [graphLevel] = applyLocalInhibition(graphLevel, options, levelItr)
    % Calculate edge radius.
    scale = (1/options.scaling)^(levelItr-1);
    neighborhood = fix(options.edgeRadius * scale);

    imageIds = [graphLevel.imageId];
    numberOfImages = max(imageIds);
    numberOfNodes = numel(graphLevel);
    discardedNodes = ones(numberOfNodes,1)>0;
    
    % Get node coordinates.
    nodeCoords = zeros(numberOfNodes,2);
    for nodeItr = 1:numberOfNodes
        nodeCoords(nodeItr,:) = graphLevel(nodeItr).position;
    end

    %% For faster calculation, we reverse array of leaf nodes.
    maxLeafNode = max([graphLevel.leafNodes]);
    leafParentArr = cell(maxLeafNode,1);
    for nodeItr = 1:numberOfNodes
        numberOfLeafNodes = numel(graphLevel(nodeItr).leafNodes);
        for leafNodeItr = 1:numberOfLeafNodes
            leafParentArr(graphLevel(nodeItr).leafNodes(leafNodeItr)) = {[ ...
                leafParentArr{(graphLevel(nodeItr).leafNodes(leafNodeItr))}, nodeItr]};
        end
    end
    
    %% Go over each node and check neighboring nodes for novelty introduced. Eliminate weak ones.
    nodeOffset = 0;
    for imageId = 1:numberOfImages
        imageNodeIdx = imageIds == imageId;
        imageGraphLevel = graphLevel(1,imageNodeIdx);
        imageNodeCoords = nodeCoords(imageNodeIdx,:);
        numberOfNodesInImage = numel(imageGraphLevel);
        for nodeItr = 1:(numberOfNodesInImage-1)
          actualNodeItr = nodeItr + nodeOffset;
          %% If nobody has erased this node before, it has a right to be in the final graph.
          if discardedNodes(actualNodeItr) == 0
              continue;
          end
          selfLeafNodes = imageGraphLevel(nodeItr).leafNodes;
          overlappingNodes = unique([leafParentArr{selfLeafNodes}]) - nodeOffset;
          overlappingNodes = overlappingNodes(overlappingNodes > nodeItr);

          %% Get each neighboring node.
          possibleAdjacentNodeCoords = imageNodeCoords(overlappingNodes,:);
          numberOfPossibleAdjacentCoords = numel(overlappingNodes);
          thisNodeCoords = imageNodeCoords(nodeItr,:);
          centerArr = repmat(thisNodeCoords, numberOfPossibleAdjacentCoords, 1);
          distances = sqrt(sum((centerArr - possibleAdjacentNodeCoords).^2, 2));
          adjacentNodes = overlappingNodes(distances <= neighborhood); 
          
          %% Go over each adjacent node, and apply inhibition if their leaf nodes are too common under current novelty threshold.
          for adjNodeItr = 1:numel(adjacentNodes)
              adjLeafNodes = imageGraphLevel(adjacentNodes(adjNodeItr)).leafNodes;
              numberOfAdjLeafNodes = numel(adjLeafNodes);
              intersectingLeafNodes = intersect(selfLeafNodes, adjLeafNodes);
              %% Check for novelty here. If new nodes are introduced, everything is good. Else, remove the adjacent node.
              if ((numberOfAdjLeafNodes-numel(intersectingLeafNodes))/numberOfAdjLeafNodes) < options.noveltyThr
                  discardedNodes(adjacentNodes(adjNodeItr) + nodeOffset) = 0;
              end
          end
        end
        nodeOffset = nodeOffset + numberOfNodesInImage;
    end
    graphLevel = graphLevel(:,discardedNodes);
end