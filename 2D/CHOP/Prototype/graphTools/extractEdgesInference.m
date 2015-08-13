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
function [nodeEdges] = extractEdgesInference(nodes, modes, leafNodeArr, firstLevelAdjNodes, options, currentLevelId)
    %% Function initializations, reading data from main graph.
    % Calculate edge radius.
    scale = (1/options.scaling)^(currentLevelId-1);
    neighborhood = floor(options.edgeRadius * scale);
    numberOfNodes = size(nodes,1);
    nodeIds = nodes(:,1);
    nodeCoords = nodes(:,2:3);
    edgeIdMatrix = options.edgeIdMatrix;
    halfMatrixSize = (options.edgeQuantize+1)/2;
    matrixSize = [options.edgeQuantize, options.edgeQuantize];
    downsampleRatio = floor((options.edgeQuantize-1)/2) / neighborhood;
    edgeType = options.edgeType;
    
    %% Program options into variables.
    if options.fastInference
        edgeNoveltyThr = 1-options.edgeNoveltyThr;
    else
        edgeNoveltyThr = 1;
    end
    maxNodeDegree = options.maxNodeDegree;
    % Eliminate low-scored adjacency links to keep the graph degree at a constant level.
    averageNodeDegree = maxNodeDegree;
    nodeEdges = cell(numberOfNodes,1);
    
    % If there are no nodes in this image, move on.
    if numberOfNodes == 0
       nodeEdges = [];
       return; 
    end

    % Get data structures containing information about the nodes in this image.
    maxSharedLeafNodes = cellfun(@(x) numel(x), leafNodeArr) * edgeNoveltyThr;
    curAdjacentNodes = cell(numberOfNodes,1);

    %% Find all edges within this image.
    for nodeItr = 1:numberOfNodes
       centerArr = repmat(nodeCoords(nodeItr,:), numberOfNodes,1);
       distances = sqrt(sum((nodeCoords - centerArr).^2, 2));
       adjacentNodes = distances <= neighborhood;

       %% Check for edge novelty.
       adjacentNodeIdx = find(adjacentNodes);
       if currentLevelId > 1
           centerLeafNodes = leafNodeArr{nodeItr};
           commonLeafCounts = cellfun(@(x) sum(ismembc(x, centerLeafNodes)), leafNodeArr(adjacentNodes));
           novelNodes = (commonLeafCounts <= maxSharedLeafNodes(adjacentNodes))';
           adjacentNodes = adjacentNodeIdx(novelNodes);
       else
           adjacentNodes = adjacentNodeIdx;
       end
       adjacentNodes = adjacentNodes(adjacentNodes~=nodeItr);
       
      %% A further check is done to ensure boundary/surface continuity.
       if strcmp(edgeType, 'continuity') && currentLevelId > 1
           centerAdjNodes = sort(cat(1, firstLevelAdjNodes{centerLeafNodes}));
           adjLeafNodes = leafNodeArr(adjacentNodes);
           validAdjacentNodes = cellfun(@(x) nnz(ismembc(x, centerAdjNodes)), adjLeafNodes) > 0;
           adjacentNodes = adjacentNodes(validAdjacentNodes);
       end     

       %% Eliminate adjacent which are far away, if the node has too many neighbors.
       % Calculate scores (distances).
       scores = distances(adjacentNodes);

       % Eliminate nodes having lower scores.
       if numel(adjacentNodes)>averageNodeDegree
            [idx] = getLargestNElements(scores, averageNodeDegree);
            adjacentNodes = adjacentNodes(idx);
       end

       %% Assign final adjacent nodes.
       curAdjacentNodes(nodeItr) = {[repmat(int32(nodeItr), numel(adjacentNodes),1), adjacentNodes]}; 
    end

    % Get rid of empty entries in curAdjacentNodes.
    nonemptyCurAdjacentNodeIdx = cellfun(@(x) ~isempty(x), curAdjacentNodes);
    curAdjacentNodes = curAdjacentNodes(nonemptyCurAdjacentNodeIdx);

    allEdges = cat(1, curAdjacentNodes{:});
    numberOfAllEdges = size(allEdges,1);

    if numberOfAllEdges == 0
       return;
    end

    node1Labels = nodeIds(allEdges(:,1));
    node2Labels = nodeIds(allEdges(:,2));
    node1Coords = nodeCoords(allEdges(:,1),:);
    node2Coords = nodeCoords(allEdges(:,2),:);
    edgeCoords = node2Coords - node1Coords;

    % Remove edges between overlapping (same position) nodes having same id.
    labelEqualityArr = node1Labels == node2Labels;
    validEdges = ~(labelEqualityArr & (edgeCoords(:,1) == 0 & edgeCoords(:,2) == 0));

    % If receptive fields are used, every edge is directed.
    directedArr = ones(nnz(validEdges),1, 'int32');

    % Update data structures based on removed edges.
    allEdges = allEdges(validEdges,:);
    edgeCoords = double(edgeCoords(validEdges,:));
    normalizedEdgeCoords = fix(fix(downsampleRatio * edgeCoords) + halfMatrixSize);

    % Double check in order not to go out of bounds.
    normalizedEdgeCoords(normalizedEdgeCoords < 1) = 1;
    normalizedEdgeCoords(normalizedEdgeCoords > matrixSize(1)) = matrixSize(1);

    %Find edge labels.
    matrixInd = sub2ind(matrixSize, normalizedEdgeCoords(:,1), normalizedEdgeCoords(:,2));
    edgeIds = edgeIdMatrix(matrixInd);
    edges = [allEdges(:,1:2), edgeIds, directedArr];

    % Due to some approximations in neighborhood calculations, some edgeIds might 
    % have been assigned as 0. Eliminate such cases.
    edges = edges(edgeIds>0,:);
    
    %% Finally, we assign modes to the edges.
     % Edges empty, do nothing.
     if ~isempty(edges)
          allEdgeCoords = options.edgeCoords;
          edgeNodeLabels = nodeIds(edges(:, 1:2));
          if size(edgeNodeLabels,2) ~= 2
               edgeNodeLabels = edgeNodeLabels';
          end

          % Assign each edge a to a relevant mode.
          validEdgeIdx = ones(size(edges,1),1) > 0;
          for edgeItr = 1:size(edges,1)
               relevantModes = modes(modes(:, 1) == edgeNodeLabels(edgeItr,1) & modes(:,2) == edgeNodeLabels(edgeItr,2), :);

               % If no modes exist, we delete edges.
               if isempty(relevantModes)
                    validEdgeIdx(edgeItr) = 0;
                    continue;
               end

               % Go through every relevant mode, and get the most likely
               % distribution's id.
               probs = zeros(size(relevantModes,1),1);
               for relevantModeItr = 1:size(relevantModes,1)
                    try
                         probs(relevantModeItr) = mvnpdf(allEdgeCoords(edges(edgeItr,3), :), relevantModes(relevantModeItr,4:5), [relevantModes(relevantModeItr,6:7); relevantModes(relevantModeItr,8:9)]);
                    catch %#ok<CTCH>
                         probs(relevantModeItr) = mvnpdf(allEdgeCoords(edges(edgeItr,3), :), relevantModes(relevantModeItr,4:5), relevantModes(relevantModeItr,[6, 9]));
                    end
               end
               [~, newLabelIdx] = max(probs);
               newLabel = int32(relevantModes(newLabelIdx(1),3));
               edges(edgeItr,3) = newLabel;
          end
          edges = edges(validEdgeIdx, :);
     end

   %% Assign all edges to their respective nodes in the final graph.
   if ~isempty(edges)
       for nodeItr = 1:numberOfNodes
           edgeIdx = edges(:,1) == (nodeItr);
           if nnz(edgeIdx) > 0
               nodeEdges(nodeItr) = {edges(edgeIdx,2:3)};
           end
       end
   end;
end
