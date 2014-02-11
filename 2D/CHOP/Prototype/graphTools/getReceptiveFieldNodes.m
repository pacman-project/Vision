%> Name: getReceptiveFieldNodes
%>
%> Description: Given a list of nodes, this function creates subsets of all
%> nodes in accordance receptive fields, i.e. candidate bounding boxes in 
%> which compositions can be searched for.
%>
%> @param nodes: Set of nodes in format: {nodeId, position, imageId, rfId, isCenter; ...]
%> The receptive field id (rfId) and center node flag (isCenter) is empty initially. 
%> @param currentLevelId: Level id of the current level (what an
%> explanation).
%> @param options: Program options.
%> 
%> @retval extendedNodes Final set of nodes. Each set's nodes has a different rfId.
%> Moreover, in each set, only one node is marked as the center node.
%> Which, obviously, is the center node itself ^^.
%> @retval receptiveFieldNodes
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 02.02.2014
function [ extendedNodes, receptiveFieldNodes] = getReceptiveFieldNodes( nodes, currentLevelId, options )
    % If no receptive field is required, do not process anything, and move
    % on.
    if ~options.useReceptiveField
       extendedNodes = nodes;
       [extendedNodes(:,4:5)] = deal({1});
       extendedNodes(:,6) = num2cell((1:size(nodes,1))');
       receptiveFieldNodes = 1:1:size(nodes,1);
       return;
    end
    
    % Initialize the size of the receptive field, along with other
    % variables.
    scale = (1/options.scaling)^(currentLevelId-1);
    neighborhood = fix(options.edgeRadius * scale);
%    receptiveFieldSize = fix(options.receptiveFieldSize * scale);
%    halfReceptiveFieldSize = floor(receptiveFieldSize/2);
    receptiveFieldNodeSets = cell(size(nodes,1),1);
    receptiveFieldIds = cell(size(nodes,1),1);
    nodeCoords = cell2mat(nodes(:,2));
    imageIds = cell2mat(nodes(:,3));
    centerNodes = zeros(size(nodes,1),1);
    
    %% From each image's (say, n) nodes, we center our receptive field around each node and select all nodes within the field.
    % This is the receptive field implementation.
    setNodeCount = 0;
    nodeOffset = 0;
    for imageItr = 1:max(imageIds)
        imageNodeIdx = find(imageIds == imageItr);
        imageCoords = nodeCoords(imageNodeIdx,:);
        numberOfImageNodes = numel(imageNodeIdx);
        
        % Select each node as the center, and get nearby nodes, which form
        % a set corresponding to a receptive field. Save receptive field
        % ids and center nodes.
        for nodeItr = 1:numberOfImageNodes
            centerArr = repmat(imageCoords(nodeItr,:), numberOfImageNodes, 1);
            distances = sqrt(sum((centerArr - imageCoords).^2,2));
            nodeSet = find(distances <= neighborhood);
%             nodeSet = find(imageCoords(:,1) >= (imageCoords(nodeItr,1)-halfReceptiveFieldSize) & ...
%                 imageCoords(:,1) <= (imageCoords(nodeItr,1)+halfReceptiveFieldSize) & ...
%                 imageCoords(:,2) >= (imageCoords(nodeItr,2)-halfReceptiveFieldSize) & ...
%                 imageCoords(:,2) <= (imageCoords(nodeItr,2)+halfReceptiveFieldSize));
            receptiveFieldNodeSets(nodeOffset + nodeItr) = {nodeSet + nodeOffset}; 
            receptiveFieldIds(nodeOffset + nodeItr) = {repmat(nodeItr, numel(nodeSet),1)};
            centerNodes(nodeOffset + nodeItr) = setNodeCount + find(nodeSet == nodeItr, 1, 'first');
            setNodeCount = setNodeCount + numel(nodeSet);
        end
        nodeOffset = nodeOffset + numberOfImageNodes;
    end
    
    %% Write list of nodes separately to the output node list.
    receptiveFieldNodes = cell2mat(receptiveFieldNodeSets);
    extendedNodes = nodes(receptiveFieldNodes,:);
    extendedNodes(:,4) = mat2cell(cell2mat(receptiveFieldIds), ones(setNodeCount, 1));
    extendedNodes(setdiff(1:setNodeCount, centerNodes),5) = {0};
    extendedNodes(centerNodes,5) = {1};
    extendedNodes(:,6) = mat2cell(receptiveFieldNodes, ones(setNodeCount, 1));
end

