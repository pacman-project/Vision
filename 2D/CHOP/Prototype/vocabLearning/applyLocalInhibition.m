%> Name: applyLocalInhibition
%>
%> Description: Applies local inhibition of nodes to eliminate those who do
%> not introduce novelty, into the graph. In this context, novelty is
%> measured by the ratio of leaf nodes not existing in the main graph over
%> all leaf nodes of the tested node. In case of an overlap, the node whose
%> label id is better wins.
%>
%> @param graphLevel The current graph level, ordered by labelIds.
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
    numberOfNodes = numel(graphLevel);
    discardedNodes = ones(numberOfNodes,1)>0;
    
    % Get node coordinates.
    nodeCoords = zeros(numberOfNodes,2);
    for nodeItr = 1:numberOfNodes
        nodeCoords(nodeItr,:) = graphLevel(nodeItr).position;
    end
    
    %% Go over each node and check neighboring nodes for novelty introduced. Eliminate weak ones.
    for nodeItr = 1:numberOfNodes
      %% If nobody has erased this node before, it has a right to be in the final graph.
      selfLeafNodes = graphLevel(nodeItr).leafNodes;
        
      %% TODO: Optimize this part to avoid redundant calculation.
      % Get each neighboring node.
      imageNodes = imageIds == imageIds(nodeItr);
      thisNodeCoords = nodeCoords(nodeItr,:);
      imageNodeCoords = nodeCoords(imageNodes,:);
      imageNodeCount = size(imageNodeCoords,1);
      centerArr = [ones(imageNodeCount,1) * thisNodeCoords(1), ...
          ones(imageNodeCount,1) * thisNodeCoords(2)];
      distances = sqrt(sum((centerArr - imageNodeCoords).^2, 2));
      adjacentNodes = find(distances <= neighborhood); 
      
      % Remove self node pair.
      adjacentNodes = adjacentNodes(adjacentNodes ~= nodeItr);
      for adjNodeItr = 1:numel(adjacentNodes)
          adjLeafNodes = graphLevel(adjacentNodes(adjNodeItr)).leafNodes;
          numberOfAdjLeafNodes = numel(adjLeafNodes);
          intersectingLeafNodes = intersect(selfLeafNodes, adjLeafNodes);
          
          %% Check for novelty here. If new nodes are introduced, everything is good. Else, remove the adjacent node.
          if ((numberOfAdjLeafNodes-numel(intersectingLeafNodes))/numberOfAdjLeafNodes) < options.noveltyThr
              discardedNodes(nodeItr) = 0;
          end
      end
    end
    graphLevel = graphLevel(:,discardedNodes);
end