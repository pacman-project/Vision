%> Name: generateLevels
%>
%> Description: Given nodes and edges, this function generates the current 
%> vocabulary and object graph levels. It only initializes the data
%> structures with few of their fields filled.
%>
%> @param nodes The node list including label id, position, image id,
%> image-wise receptive field id and center flag for each node. If no
%> receptive fields are used, center flag is true for each node.
%> @param options Program options.
%>
%> @retval vocabLevel The vocabulary level.
%> @retval graphLevel The object graph level.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 03.02.2014
function [vocabLevel, graphLevel] = generateLevels(nodes, options)
    % Allocate space for both levels.
    vocabLevel(options.numberOfFilters) = options.vocabNode;
    graphLevel(size(nodes,1)) = options.graphNode;
    
    % Fill vocabulary level with labels.
    labelArr = num2cell(1:6);
    [vocabLevel.label] = deal(labelArr{:});
    
    % Fill object graph level.
    [graphLevel.labelId, graphLevel.position, graphLevel.imageId] = deal(nodes{:});
    
    % Set the sign of all nodes to 1. When negative graphs are introduced,
    % this part should CHANGE.
    [graphLevel.sign] = deal(true);
    
    % Add leaf nodes.
    leafNodes = num2cell(1:size(nodes,1));
    [graphLevel.leafNodes] = deal(leafNodes{:});
end