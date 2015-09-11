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
%> Ver 1.1 on 01.09.2014 Adding negative signs.
function [vocabLevel, graphLevel] = generateLevels(nodes, activations, signs, options)
    % Allocate space for both levels.
    vocabLevel(options.numberOfFilters) = VocabNode();
    graphLevel(size(nodes,1)) = GraphNode();
    
    % Fill vocabulary level with labels.
    labelArr = num2cell(int32(1:options.numberOfFilters));
    [vocabLevel.label] = deal(labelArr{:});
    [vocabLevel.numberOfProbabilisticChoices] = deal(int32(1));
    clear labelArr;
    
    % Fill object graph level.
    assignedNodes = mat2cell(nodes, ones(size(nodes,1),1), [1 2 1]);
    [graphLevel.labelId, graphLevel.position, graphLevel.imageId] = deal(assignedNodes{:});
    assignedPrecisePos = mat2cell(single(nodes(:,2:3)), ones(size(nodes,1),1), 2);
    [graphLevel.precisePosition] = deal(assignedPrecisePos{:});
    clear assignedNodes;
    
    % Set the signs.
    signArr = num2cell(signs);
    [graphLevel.sign] = deal(signArr{:});
    clear signArr;
    
    % Set activation values.
    activationArr = num2cell(activations);
    [graphLevel.activation] = deal(activationArr{:});
    
    % Add leaf nodes.
    leafNodes = num2cell(int32(1:size(nodes,1)));
    [graphLevel.leafNodes] = deal(leafNodes{:});
    clear leafNodes;
end