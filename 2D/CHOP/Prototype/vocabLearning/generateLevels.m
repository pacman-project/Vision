%> Name: generateLevels
%>
%> Description: Given nodes and edges, this function generates the current 
%> vocabulary and object graph levels. It only initializes the data
%> structures with few of their fields filled.
%>
%> @param nodes The node list including label id, position, image id,
%> image-wise receptive field id and center flag for each node. If no
%> receptive fields are used, center flag is true for each node.
%> @param leafNodes The leaf node list including label id, position, image
%> id of each leaf node. Also referred as level 0 nodes.
%> @param options Program options.
%>
%> @retval vocabLevel The vocabulary level.
%> @retval graphLevel The object graph level.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 03.02.2014
function [vocabLevel, graphLevel] = generateLevels(nodes, leafNodes, options)
    % Allocate space for both levels.
    vocabLevel(options.numberOfFilters) = options.vocabNode;
    graphLevel(size(nodes,1)) = options.graphNode;
    
    % Fill vocabulary level with labels.
    labelArr = num2cell(1:6);
    [vocabLevel.label] = deal(labelArr{:});
    
    % Fill object graph level.
    [graphLevel.labelId, graphLevel.position, graphLevel.imageId, ...
        graphLevel.rfId, graphLevel.isCenter, graphLevel.realNodeId] = deal(nodes{:});
    
    % Set the sign of all nodes to 1. When negative graphs are introduced,
    % this part should CHANGE.
    [graphLevel.sign] = deal(1);
    
    leafNodePositions = cell2mat(leafNodes(:,2));
    numberOfLeafNodes = size(leafNodePositions,1);
    leafNodeImageIds = cell2mat(leafNodes(:,3));
    % Add leaf nodes and edge info.
    for instanceItr = 1:size(nodes,1)
       leafNode = find(sum(abs(leafNodePositions - ...
        repmat(graphLevel(instanceItr).position, numberOfLeafNodes,1)),2)==0 & ...
        leafNodeImageIds == graphLevel(instanceItr).imageId, 1, 'first');
    
        graphLevel(instanceItr).leafNodes = leafNode;
%        
%        edges = ismember(edges(:,1:2), instanceItr);
%        % get non-zero rows of edges
%        selfEdges = edges(edges(:,1) | edges(:,2),:);
%        if numel(selfEdges)>0
%             graphLevel(instanceItr).adjInfo = selfEdges;
%        end
    end
end