%> Name: createEdgesWithLabels
%>
%> Description: Given the node list, the aim of this
%> function is to detect spatial distributions by analyzing 2-D spatial
%> arrangements of node types. For example, if there are 3 distinct
%> configurations of a node type pair, the edges in between are represented 
%> with 3 categories.
%>
%> @param mainGraph The object graphs' data structure.
%> @param options Program options.
%> @param currentLevel The current scene graph level.
%> @param modes If empty, new modes are to be learned. If not, edge labels
%>      will be formed depending on existing modes.
%> @param hMatrix histogram matrix helping to label edges in 'hist' mode.
%> 
%> @retval edges Edges are of the form: [ node1, node2, mode, directed;
%>                                        node1, node2, mode, directed;
%>                                      ...]
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 04.12.2013
%> Separate mode learning from this function on 28.01.2014
function [mainGraph] = createEdgesWithLabels(mainGraph, options, currentLevelId, modes, hMatrix)
    %% Function initializations, reading data from main graph.
    % Calculate edge radius.
    scale = (1/options.scaling)^(currentLevelId-1);
    neighborhood = fix(options.edgeRadius * scale);
    
    currentLevel = mainGraph{currentLevelId};
    nodeIds = [currentLevel.labelId]';
    nodeCoords = cat(1, currentLevel.position);
    imageIds = [currentLevel.imageId]';
    
    %% Program options into variables.
    edgeNoveltyThr = 1-options.edgeNoveltyThr;
    maxNodeDegree = options.maxNodeDegree;
    % Eliminate low-scored adjacency links to keep the graph degree at a constant level.
    averageNodeDegree = maxNodeDegree;
    
    property = options.property;
    %% Calculate size of the histogram matrix.
    sizeHMatrix = size(hMatrix,1);
    halfSizeHMatrix = floor(sizeHMatrix/2);
    
    %% Get relevant pair-wise distributions (nodes)
    if ~isempty(modes) && numel(modes) >= currentLevelId
        currentModes = modes{currentLevelId};
        currentModesFirstNodes = currentModes(:,1);
        currentModesSecNodes = currentModes(:,2);
        currentModesPos = currentModes(:,3:4);
    else
        currentModesFirstNodes = [];
        currentModesSecNodes = [];
        currentModesPos = [];
    end
    
    %% Put each image's node set into a different bin.
    numberOfImages = max(imageIds);
    imageGraphNodeSets = cell(numberOfImages, 1);
    imageNodeIdArr = cell(numberOfImages,1);
    imageNodeCoordArr = cell(numberOfImages,1);
    for imageItr = 1:max(imageIds)
       imageNodeIdx = imageIds == imageItr;
       imageGraphNodeSets(imageItr) = {currentLevel(imageIds == imageItr)};
       imageNodeIdArr(imageItr) = {nodeIds(imageNodeIdx)};
       imageNodeCoordArr(imageItr) = {nodeCoords(imageNodeIdx,:)};
    end
    
    
    %% Process each image separately (and in parallel)
    parfor imageItr = 1:numberOfImages
        imageNodeOffset = numel(find(imageIds < imageItr));
        imageNodeIdx = imageIds == imageItr;

        % If there are no nodes in this image, move on.
        if nnz(imageNodeIdx) == 0
           continue; 
        end
        
        % Get data structures containing information about the nodes in this image.
        curNodeIds = imageNodeIdArr{imageItr};
        curNodeCoords = imageNodeCoordArr{imageItr};
        imageNodeIdx = find(imageNodeIdx)';
        numberOfNodes = numel(imageNodeIdx);
        curGraphNodes = imageGraphNodeSets{imageItr};
        curLeafNodes = {curGraphNodes.leafNodes};
        maxSharedLeafNodes = cellfun(@(x) numel(x), curLeafNodes) * edgeNoveltyThr;
        curAdjacentNodes = cell(numberOfNodes,1);
        
        %% Find all edges within this image.
        for nodeItr = 1:numberOfNodes
           centerArr = repmat(curNodeCoords(nodeItr,:), numberOfNodes,1);
           distances = sqrt(sum((centerArr - curNodeCoords).^2, 2));
           adjacentNodes = distances <= neighborhood;
           
           %% Check for edge novelty.
           adjacentNodeIdx = find(adjacentNodes);
           if currentLevelId > 1
               centerLeafNodes = curLeafNodes{nodeItr};
               commonLeafCounts = cellfun(@(x) sum(ismembc(x, centerLeafNodes)), curLeafNodes(adjacentNodes));
               novelNodes = (commonLeafCounts <= maxSharedLeafNodes(adjacentNodes))';
               adjacentNodes = adjacentNodeIdx(novelNodes);
           else
               adjacentNodes = adjacentNodeIdx;
           end
           adjacentNodes = adjacentNodes(adjacentNodes~=nodeItr);
           
           %% Eliminate adjacent which are far away, if the node has too many neighbors.
           % Calculate scores (distances).
           scores = distances(adjacentNodes);
           
           % Eliminate nodes having lower scores.
           if numel(adjacentNodes)>averageNodeDegree
                [idx] = getSmallestNElements(scores, averageNodeDegree);
                adjacentNodes = adjacentNodes(idx);
           end
           
           %% Assign final adjacent nodes.
           curAdjacentNodes(nodeItr) = {[repmat(nodeItr, numel(adjacentNodes),1), adjacentNodes]}; 
        end
        
        % Get rid of empty entries in curAdjacentNodes.
        nonemptyCurAdjacentNodeIdx = cellfun(@(x) ~isempty(x), curAdjacentNodes);
        curAdjacentNodes = curAdjacentNodes(nonemptyCurAdjacentNodeIdx);
        
        allEdges = cat(1, curAdjacentNodes{:});
        numberOfAllEdges = size(allEdges,1);
        
        if numberOfAllEdges == 0
           continue;
        end
        
        node1Labels = curNodeIds(allEdges(:,1));
        node2Labels = curNodeIds(allEdges(:,2));
        node1Coords = curNodeCoords(allEdges(:,1),:);
        node2Coords = curNodeCoords(allEdges(:,2),:);
        edgeCoords = node1Coords - node2Coords;
        
        % Remove edges between overlapping (same position) nodes having same id.
        labelEqualityArr = node1Labels == node2Labels;
        validEdges = ~(labelEqualityArr & (edgeCoords(:,1) == 0 & edgeCoords(:,2) == 0));

        % If receptive fields are used, every edge is directed.
        directedArr = ones(nnz(validEdges),1);
        
        % Update data structures based on removed edges.
        allEdges = allEdges(validEdges,:);
        edgeCoords = edgeCoords(validEdges,:);
        numberOfAllEdges = nnz(validEdges);
        node1Labels = node1Labels(validEdges,:);
        node2Labels = node2Labels(validEdges,:);
        
        %% Based on the edge type, assign edge labels here.
        if strcmp(property, 'hist')
            normalizedEdgeCoords = round((edgeCoords / neighborhood)*(halfSizeHMatrix-2)) + halfSizeHMatrix;
            hMatrixInd = sub2ind([sizeHMatrix, sizeHMatrix], normalizedEdgeCoords(:,1), normalizedEdgeCoords(:,2));
            edgeIds = hMatrix(hMatrixInd);
        elseif strcmp(property, 'mode');
            %% Process each node in this image. 
           edgeIds = zeros(numberOfAllEdges,1);
           for edgeItr = 1:numberOfAllEdges
                % Estimate mode of this edge.
                applicableModeIdx = find(currentModesFirstNodes == node1Labels(edgeItr) & ...
                currentModesSecNodes == node2Labels(edgeItr));
                applicableModes = currentModesPos(applicableModeIdx,:);
                centerArr2 = repmat(edgeCoords(edgeItr,:), size(applicableModes,1),1);
                distances = sqrt(sum((centerArr2 - applicableModes).^2, 2));
                [~, minDist] = min(distances);

                % If a valid mode exists, assign its label.
                if numel(minDist)>0
                    edgeIds(edgeItr) = applicableModeIdx(minDist(1));
                else
                    edgeIds(edgeItr) = 0;
                end
           end
        else
            edgeIds = ones(numberOfAllEdges,1);
        end
        
        edges = [allEdges(:,1:2) + imageNodeOffset, edgeIds, directedArr];
        
        % Due to some approximations in mode calculations, some edgeIds might 
        % have been assigned as 0. Eliminate such cases.
        edges = edges(edgeIds>0,:);
        
       %% Assign all edges to their respective nodes in the final graph.
       
       if ~isempty(edges)
           for nodeItr = 1:numberOfNodes
               curGraphNodes(nodeItr).adjInfo = edges(edges(:,1) == (nodeItr + imageNodeOffset),:);
           end
       end
       imageGraphNodeSets(imageItr) = {curGraphNodes};
    end
    currentLevel = [imageGraphNodeSets{:}];
    mainGraph(currentLevelId) = {currentLevel};
end
