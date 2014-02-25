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
function [mainGraph] = createEdgesWithLabels(mainGraph, options, currentLevelId, modes)
    %% Function initializations, reading data from main graph.
    % Calculate edge radius.
    scale = (1/options.scaling)^(currentLevelId-1);
    neighborhood = fix(options.edgeRadius * scale);
    
    currentLevel = mainGraph{currentLevelId};
    nodeIds = [currentLevel.labelId]';
    nodeCoords = cat(1, currentLevel.position);
    imageIds = [currentLevel.imageId]';
    
    %% Program options into variables.
    edgeNoveltyThr = options.edgeNoveltyThr;
    useReceptiveField = options.useReceptiveField;
    maxNodeDegreeLevel1 = options.maxNodeDegreeLevel1;
    maxNodeDegree = options.maxNodeDegree;
    property = options.property;
    %% Read histogram matrix to be used in 'hist' type edge calculations.
%    if strcmp(options.property, 'hist')
    load('hMatrix.mat', 'hMatrix'); 
    sizeHMatrix = size(hMatrix,1);
    halfSizeHMatrix = floor(sizeHMatrix/2);
%    end
    
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
    imageGraphNodeSets = cell(max(imageIds), 1);
    for imageItr = 1:max(imageIds)
       imageGraphNodeSets(imageItr) = {currentLevel(imageIds == imageItr)};
    end
    
    %% Process each image separately (and in parallel)
    parfor imageItr = 1:max(imageIds)
        imageNodeIdx = imageIds == imageItr;
        nodeOffset = numel(find(imageIds < imageItr));

        % If there are no nodes in this image, move on.
        if nnz(imageNodeIdx) == 0
           continue; 
        end
        
        % Get data structures corresponding to nodes in this image.
        curNodeIds = nodeIds(imageNodeIdx);
        curNodeCoords = nodeCoords(imageNodeIdx,:);
        imageNodeIdx = find(imageNodeIdx)';
        numberOfNodes = numel(imageNodeIdx);
        curGraphNodes = imageGraphNodeSets{imageItr};
        curLeafNodes = {curGraphNodes.leafNodes};

        %% Process each node in this image. 
        for nodeItr = 1:numel(imageNodeIdx)
           %% Depending on the edge type, find a list of nodes which are neighbors of this node.
           centerArr = repmat(curNodeCoords(nodeItr,:), numberOfNodes,1);
           distances = sqrt(sum((centerArr - curNodeCoords).^2, 2));
           nearbyNodes = distances <= neighborhood;

           % Check adjacent nodes to see if they are among nodes
           % introducing enough novelty with respect to each other. 
           % options.edgeNoveltyThr provides a lower bound for the
           % novelty each edge should introduce.
           if currentLevelId > 1
               centerLeafNodes = curGraphNodes(nodeItr).leafNodes;
               leafCounts = cellfun(@(x) numel(x), curLeafNodes);
               diffLeafCounts = cellfun(@(x) numel(ismembc(x, centerLeafNodes)), curLeafNodes);
               novelNodes = (diffLeafCounts >= edgeNoveltyThr*leafCounts)';
               adjacentNodes = novelNodes & nearbyNodes;
           else
               adjacentNodes = nearbyNodes;
           end

           % 'distance' type score calculation
           adjacentNodes = find(adjacentNodes);

           % Looking for closest nodes to connect!
           scores = distances(adjacentNodes); 

           % Let's remove center node from its neighbors.
           otherAdjacentNodes = adjacentNodes ~= nodeItr;
           adjacentNodes = adjacentNodes(otherAdjacentNodes);
           scores = scores(otherAdjacentNodes);

           % If no adjacent nodes exist, exit.
           if isempty(adjacentNodes)
               continue;
           end

           %% If using receptive fields, there is no limit on the number of neighbors.
           if ~useReceptiveField
               % Get valid adjacent nodes with larger node ids.
               validAdjacentNodes = curNodeIds(adjacentNodes) >= curNodeIds(nodeItr);
               adjacentNodes = adjacentNodes(validAdjacentNodes);
               scores = scores(validAdjacentNodes);

               % Eliminate low-scored adjacency links to keep the graph degree at a constant level.
               if currentLevelId == 1
                   averageNodeDegree = maxNodeDegreeLevel1;
               else
                   averageNodeDegree = maxNodeDegree;
               end
               if numel(adjacentNodes)>averageNodeDegree
                    [idx] = getSmallestNElements(scores, averageNodeDegree);
                    adjacentNodes = adjacentNodes(idx);
               end
           end

           if useReceptiveField
               directedArr = ones(numel(adjacentNodes),1);
           else
               samples = curNodeCoords(adjacentNodes,:) - repmat(curNodeCoords(nodeItr,:), numel(adjacentNodes),1);
               sampleSums = sum(samples,2);
               
               % Check which edges are on the right side of separating line.
               validSamples = ~(curNodeIds(adjacentNodes) == curNodeIds(nodeItr)) | ...
                   ((samples(:,1) >= 0 & sampleSums >= 0) | (samples(:,1) < 0 & sampleSums > 0));
               
               % If both are at the origin (same spot), this causes
               % into duplicate links in our final graph. To prevent
               % this, we only link one with smaller index to the
               % other one with larger index.
               validSamples = ~(curNodeIds(adjacentNodes) == curNodeIds(nodeItr)) | ...
                   (validSamples & (~(samples(:,1) == 0 & samples(:,2) == 0) | adjacentNodes>nodeItr));
               adjacentNodes = adjacentNodes(validSamples);
               
               % If both have same label, edge is directed. Otherwise, it
               % is undirected.
               directedArr = curNodeIds(adjacentNodes) == curNodeIds(nodeItr);
           end
           
           samples = curNodeCoords(adjacentNodes,:) - repmat(curNodeCoords(nodeItr,:), numel(adjacentNodes),1);
           %% Create an edge for each center-adjacent node pair.
           numberOfAdjacentNodes = numel(adjacentNodes);
           if strcmp(property, 'hist')
               % Assign the edges their labels based on the histogram info
               normalizedSamples = round((double(samples) / neighborhood)*halfSizeHMatrix) + halfSizeHMatrix;
               normalizedSamples(normalizedSamples < 1) = 1;
               normalizedSamples(normalizedSamples > sizeHMatrix) = sizeHMatrix;
               hMatrixInd = sub2ind([sizeHMatrix, sizeHMatrix], normalizedSamples(:,1), normalizedSamples(:,2));
               edgeIds = hMatrix(hMatrixInd);
           elseif strcmp(property, 'co-occurence')
               edgeIds = ones(numberOfAdjacentNodes,1);
           else
               edgeIds = zeros(numberOfAdjacentNodes,1);
               % Mode-based edge label assignment.
               for adjNodeItr = 1:numberOfAdjacentNodes
                    % Estimate mode of this edge.
                    applicableModeIdx = find(currentModesFirstNodes == curNodeIds(nodeItr) & ...
                    currentModesSecNodes == curNodeIds(adjacentNodes(adjNodeItr)));
                    applicableModes = currentModesPos(applicableModeIdx,:);
                    centerArr2 = repmat(samples(adjNodeItr,:), size(applicableModes,1),1);
                    distances = sqrt(sum((centerArr2 - applicableModes).^2, 2));
                    [~, minDist] = min(distances);

                    % If a valid mode exists, assign its label.
                    if numel(minDist)>0
                        edgeIds(adjNodeItr) = applicableModeIdx(minDist(1));
                    else
                        edgeIds(adjNodeItr) = 0;
                    end
               end
           end
           edges = [repmat(imageNodeIdx(nodeItr), numberOfAdjacentNodes,1), ...
               adjacentNodes + nodeOffset, ...
               edgeIds, ...
               directedArr];
           
           %% Assign all edges to the node in the given graph.
           if ~isempty(edges)
               curGraphNodes(nodeItr).adjInfo = edges;
           end
        end
        imageGraphNodeSets(imageItr) = {curGraphNodes};
    end
    currentLevel = [imageGraphNodeSets{:}];
    mainGraph(currentLevelId) = {currentLevel};
end
