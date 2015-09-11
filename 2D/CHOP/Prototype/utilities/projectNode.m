%> Name: projectNode
%>
%> Description: Given the nodes in 'nodes' array, and vocabulary descriptions 
%> in vocabulary, we backproject the nodes to image plane by recursively
%> projecting sub-nodes.
%>
%> @param nodes exported nodes in the format of 
%> [labelId, centerPosX, centerPosY, levelId]
%>   
%> @param vocabulary The vocabulary.
%> @param inhibitionRadius The radius in which we will suppress other responses.
%>  
%> @retval nodes Projected level 1 nodes. 
%>
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 12.08.2015
function [ nodes ] = projectNode( nodes, vocabulary, inhibitionRadius )
    levelItr = nodes(1,4);
    nodes = single(nodes);
    inhibitionRadius = inhibitionRadius - 1;
    inhibitionRadiusSq = inhibitionRadius^2;
    
    %% First, we recursively backproject the nodes. 
    while levelItr > 1.001
        vocabLevel = vocabulary{levelItr};
         
         % Get real ids of the parts.
         vocabLevelLabels = [vocabLevel.label];
         realIdx = zeros(max(vocabLevelLabels),1);
         for itr = 1:size(realIdx,1)
              realIdx(itr) = find(vocabLevelLabels == itr, 1, 'first');
         end
         
        newNodes = cell(size(nodes,1),1);
        for nodeItr = 1:size(nodes,1)
            vocabNode = vocabLevel(realIdx(nodes(nodeItr,1)));
            newNodeSet = zeros(numel(vocabNode.children), 4, 'single');
            newNodeSet(:, 1) = single((vocabNode.children)');
            newNodeSet(:, 4) = levelItr-1;
            newNodeSet(:, 2:3) = (vocabNode.childrenPosMean) + repmat(nodes(nodeItr, 2:3), size(newNodeSet,1), 1);
            newNodes{nodeItr} = newNodeSet;
        end
        nodes = cat(1, newNodes{:});
        levelItr = levelItr - 1;
    end
    nodes = int32(round(nodes(:,1:3)));
    
    %% Then, we perform a simple inhibition process to remove overlapping level 1 instances.
    validNodes = ones(size(nodes,1),1) > 0;
%    numberOfNodes = size(nodes,1);
%     for nodeItr = 1:(numberOfNodes-1)
%         if ~validNodes(nodeItr) 
%            continue; 
%         end
%         remainingNodes = nodes((nodeItr+1):end, :);
%         validRemainingNodes = sum((remainingNodes(:,2:3) - repmat(nodes(nodeItr,2:3), ...
%             numberOfNodes - nodeItr, 1)).^2, 2) >= inhibitionRadiusSq;
%         validNodes((nodeItr+1):end) = validNodes((nodeItr+1):end) & validRemainingNodes;
%     end
    
    nodes = nodes(validNodes, :);
end

