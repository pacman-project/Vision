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
function [nodeEdges, edgeProbs] = extractEdgesInference(nodes, modes, modeProbArr, leafNodeArr, firstLevelAdjNodes, options, currentLevelId)
    %% Function initializations, reading data from main graph.
    % Calculate edge radius.
    numberOfNodes = size(nodes,1);
    nodeIds = nodes(:,1);
    nodeCoords = nodes(:,2:3);
    edgeIdMatrix = options.edgeIdMatrix;
    halfMatrixSize = (options.receptiveFieldSize+1)/2;
    matrixSize = [options.receptiveFieldSize, options.receptiveFieldSize];
    edgeType = options.edgeType;
    edgeCoords = options.edgeCoords + halfMatrixSize;
    
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
    edgeProbs = cell(numberOfNodes,1);
    
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
       adjacentNodes = nodeCoords(:,1) > (centerArr(:,1) - halfMatrixSize) & ...
           nodeCoords(:,1) < (centerArr(:,1) + halfMatrixSize) & ...
           nodeCoords(:,2) > (centerArr(:,2) - halfMatrixSize) & ...
           nodeCoords(:,2) < (centerArr(:,2) + halfMatrixSize);

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

    % Obtain edges and count them.
    allEdges = cat(1, curAdjacentNodes{:});
    numberOfAllEdges = size(allEdges,1);

    if numberOfAllEdges == 0
       return;
    end

    node1Labels = nodeIds(allEdges(:,1));
    node2Labels = nodeIds(allEdges(:,2));
    node1Coords = nodeCoords(allEdges(:,1),:);
    node2Coords = nodeCoords(allEdges(:,2),:);
    allEdgeCoords = node2Coords - node1Coords;

    % Remove edges between overlapping (same position) nodes having same id.
    labelEqualityArr = node1Labels == node2Labels;
    validEdges = ~(labelEqualityArr & (allEdgeCoords(:,1) == 0 & allEdgeCoords(:,2) == 0));

    % If receptive fields are used, every edge is directed.
    directedArr = ones(nnz(validEdges),1, 'int32');

    % Update data structures based on removed edges.
    allEdges = allEdges(validEdges,:);
    allEdgeCoords = allEdgeCoords(validEdges,:);
    node1Labels = node1Labels(validEdges,:);
    
    % Double check in order not to go out of bounds.
    node2Labels = node2Labels(validEdges,:);
    allEdgeCoords = allEdgeCoords + halfMatrixSize;
    
    %Find edge labels.
    matrixInd = sub2ind(matrixSize, allEdgeCoords(:,1), allEdgeCoords(:,2));
    edgeIds = edgeIdMatrix(matrixInd);
    edges = [allEdges(:,1:2), edgeIds, directedArr];
    
    %% Finally, we assign modes to the edges.
     % Edges empty, do nothing.
     probArr = zeros(size(edges,1),1, 'single');
     if ~isempty(edges)
          % Assign each edge a to a relevant mode.
          for edgeItr = 1:size(edges,1)
               relevantIdx = modes(:, 1) == node1Labels(edgeItr) & modes(:,2) == node2Labels(edgeItr);
               relevantModes = modes(relevantIdx,:);
               if isempty(relevantModes)
                    continue;
               end
               relevantCoords = edgeCoords(edges(edgeItr,3),:);
               clusterProbs = modeProbArr(relevantIdx,relevantCoords(1),relevantCoords(2));
               [probability, clusterId] = max(clusterProbs);
               probArr(edgeItr) = probability;
               newLabel = int32(relevantModes(clusterId,3));
               edges(edgeItr,3) = newLabel;
          end
     end

   %% Assign all edges to their respective nodes in the final graph.
   if ~isempty(edges)
       for nodeItr = 1:numberOfNodes
           edgeIdx = edges(:,1) == (nodeItr);
           if nnz(edgeIdx) > 0
               nodeEdges(nodeItr) = {edges(edgeIdx,2:3)};
               edgeProbs(nodeItr) = {probArr(edgeIdx)};
           end
       end
   end;
end
