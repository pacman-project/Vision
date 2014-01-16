%> Name: calculateAndAssignModes
%>
%> Description: Given the node list, the aim of this
%> function is to detect spatial distributions by analyzing 2-D spatial
%> arrangements of node types. For example, if there are 3 distinct
%> configurations of a node type pair, the edges in between are represented 
%> with 3 categories.
%>
%> @param nodes The node list including label ids, positions and image ids
%> of each node.
%> @param mainGraph The object graphs' data structure.
%> @param options Program options.
%> @param currentLevel The current scene graph level.
%> @param datasetName Name of the dataset.
%> @param modes If empty, new modes are to be learned. If not, edge labels
%> will be formed depending on existing modes.
%> 
%> @retval modes The mode list representing edge categories.
%> @retval edges Edges are of the form: [ node1, node2, mode, directed;
%>                                        node1, node2, mode, directed;
%>                                      ...]
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 04.12.2013
function [modes, edges, leafNodeAdjArr] = calculateAndAssignModes(nodes, mainGraph, leafNodeAdjArr, options, currentLevel, datasetName, modes)

    % Calculate edge radius.
    scale = (1/options.scaling)^(currentLevel-1);
    neighborhood = fix(options.edgeRadius * scale);
    
    nodeIds = cell2mat(nodes(:,1));
    nodeCoords = cell2mat(nodes(:,2));
    imageIds = cell2mat(nodes(:,3));
    
    %% Get current graph level
    if strcmp(options.edgeType, 'contour') && ~isempty(mainGraph)
        currentGraphLevel = mainGraph{currentLevel};
    end
    edges = [];
    
    %% In level 1, we create first level adjacency information between nodes.
    if currentLevel==1 && strcmp(options.edgeType, 'contour')
        leafNodeAdjArr = cell(size(nodes,1),1);
    end
    
    %% Learn pair-wise node distributions (modes)
    if isempty(modes) 
        [currentModes] = learnModes(nodes, mainGraph, leafNodeAdjArr, options, currentLevel, datasetName, modes);    
        modes = currentModes;
    else
        currentModes = modes{currentLevel};
    end
    
    %% Assuming the distributions have already been learned, create edges and assign edge labels accordingly.
    nodeOffset = 0;
    for imageItr = 1:max(imageIds)
        imageNodes = nodes(imageIds == imageItr,:);
        imageCoords = nodeCoords(imageIds == imageItr,:);
        imageNodeIds = nodeIds(imageIds == imageItr,:);
        imageNodeIdx = find(imageIds == imageItr);

        % If there are no nodes in this image, move on.
        if isempty(imageNodes)
           continue; 
        end

        numberOfNodes = size(imageCoords,1);
        for nodeItr = 1:numberOfNodes
           centerArr = [ones(size(imageNodes,1),1) * imageCoords(nodeItr,1), ...
               ones(size(imageNodes,1),1) * imageCoords(nodeItr,2)];

           if strcmp(options.edgeType, 'contour') && ~isempty(mainGraph)
               % In contour type edges, we work on the first
               % level to reveal existing neighboring
               % information. High level compositions are
               % linked with their adjacency information among
               % their leaf nodes at the first level.
               adjacentNodes = zeros(numel(imageNodeIdx),1);
               selfLeafNodes = currentGraphLevel(imageNodeIdx(nodeItr)).leafNodes;
               selfLeafNeighbors = cell2mat(leafNodeAdjArr(selfLeafNodes));
               scores = zeros(numel(imageNodeIdx),1);
               for secNodeItr = 1:numel(imageNodeIdx)
                  %% Score calculation for each edge takes place in here.
                  bondStrength = sum(ismember(currentGraphLevel(imageNodeIdx(secNodeItr)).leafNodes, selfLeafNeighbors));
                  maxOverlap = numel(currentGraphLevel(imageNodeIdx(secNodeItr)).leafNodes) * (1-options.edgeNoveltyThr);
                  numberOfOverlapLeaves = sum(ismember(selfLeafNodes, currentGraphLevel(imageNodeIdx(secNodeItr)).leafNodes));
                  if currentLevel == 1 || ((bondStrength > 0 && numberOfOverlapLeaves <= maxOverlap))
                     scores(secNodeItr) = numberOfOverlapLeaves - bondStrength;
                     adjacentNodes(secNodeItr) = 1;
                  end
               end
               adjacentNodes = find(adjacentNodes);
           else
                distances = sqrt(sum((centerArr - imageCoords).^2, 2));
                adjacentNodes = find(distances <= neighborhood);
                scores = distances(adjacentNodes);  
           end
           
           otherAdjacentNodes = adjacentNodes ~= nodeItr;
           adjacentNodes = adjacentNodes(otherAdjacentNodes);
           scores = scores(otherAdjacentNodes);
           validAdjacentNodes = imageNodeIds(adjacentNodes) >= imageNodeIds(nodeItr);
           adjacentNodes = adjacentNodes(validAdjacentNodes);
           scores = scores(validAdjacentNodes);

           %% Eliminate low-scored adjacency links to keep the graph degree at a constant level.
           if currentLevel == 1
               averageNodeDegree = options.maxNodeDegreeLevel1;
           else
               averageNodeDegree = options.maxNodeDegree;
           end
           if numel(adjacentNodes)>averageNodeDegree
                [idx] = getSmallestNElements(scores, averageNodeDegree);
                adjacentNodes = adjacentNodes(idx);
           end
           
           % For level 1, add adjacency information to the
           % leafNodeAdjArr array. Please note that this
           % adjacency information is directed (a->b =/= b->a)
           if strcmp(options.edgeType, 'contour') && currentLevel == 1
               leafNodeAdjArr(imageNodeIdx(nodeItr)) = {[leafNodeAdjArr{imageNodeIdx(nodeItr)}; imageNodeIdx(adjacentNodes)]};
           end

           for adjItr = 1:numel(adjacentNodes)
               sample = imageCoords(adjacentNodes(adjItr),:) - imageCoords(nodeItr,:);
               if imageNodeIds(nodeItr) == imageNodeIds(adjacentNodes(adjItr))
                   isDirected = 1;
                   % Check if this edge is on the right side of separating line.
                   if sample(1) + sample(2) < 0
                       continue;
                   end
               else
                   isDirected = 0;
               end

               % Learn mode by assigning a cluster center to its
               % relative positioning.
               applicableModeIdx = find(currentModes(:,1) == imageNodeIds(nodeItr) & ...
                   currentModes(:,2) == imageNodeIds(adjacentNodes(adjItr)));
               applicableModes = currentModes(applicableModeIdx,3:4);
               centerArr2 = [ones(size(applicableModes,1),1) * sample(1), ones(size(applicableModes,1),1) * sample(2)];
               distances = sqrt(sum((centerArr2 - applicableModes).^2, 2));
               [~, minDist] = min(distances);
               if numel(minDist)>0
                   modeId = applicableModeIdx(minDist(1));
                   edges = [edges; [nodeItr + nodeOffset, adjacentNodes(adjItr) + nodeOffset, modeId, isDirected]];
               end
           end
        end
        nodeOffset = nodeOffset + numberOfNodes;
    end
    %% Remove duplicate neighbors from leaf node adjacency list.
    if currentLevel==1 && strcmp(options.edgeType, 'contour')
        leafNodeAdjArr = cellfun(@(x) unique(x), leafNodeAdjArr, 'UniformOutput', false);
    end
end

function [idx] = getSmallestNElements(vect, n)
    [~, sortedIdx] = sort(vect);
    idx = sortedIdx(1:n);
    restVect = vect(sortedIdx((n+1):end));
    % Add nodes which qualify, since same-valued nodes are already in idx.
    % Number of such nodes cannot exceed n.
    restOfNodes = (find(restVect == vect(idx(end)))+n);
    if numel(restOfNodes)>n
       restOfNodes = restOfNodes(1:n); 
    end
    idx = [idx; restOfNodes];
end

function [modes] = learnModes(nodes, mainGraph, leafNodeAdjArr, options, currentLevel, datasetName, modes)
    % Prevent empty cluster warnings in kmeans.
    w = warning('off', 'all');
    
    % Calculate edge radius.
    scale = (1/options.scaling)^(currentLevel-1);
    neighborhood = fix(options.edgeRadius * scale);
    
    % Set initial data structures for processing 
    nodeIds = cell2mat(nodes(:,1));
    nodeCoords = cell2mat(nodes(:,2));
    imageIds = cell2mat(nodes(:,3));
    numberOfNodes = max(nodeIds);
    if strcmp(options.edgeType, 'contour') && ~isempty(mainGraph)
        currentGraphLevel = mainGraph{currentLevel};
    end
    
    %% For each node pair, get the 2-D distribution of samples.
    % Right now, we only work with 2-dimensional spatial relations.

    modeOffset = 0;

    for node1 = 1:numberOfNodes
       for node2 = node1:numberOfNodes
           % Get first type of nodes.
           firstNodeIdx = find(nodeIds==node1);
           firstNodeCoords = nodeCoords(firstNodeIdx, :);
           firstNodeImageIds = imageIds(firstNodeIdx, :);

           % Get second type of nodes.
           secondNodeIdx = find(nodeIds==node2);
           secondNodeCoords = nodeCoords(secondNodeIdx, :);
           secondNodeImageIds = imageIds(secondNodeIdx, :);

           % Allocate space for pair-wise distribution data
           samples = [];
           % Find all edges in between.
           for imageId = 1:max(firstNodeImageIds)
               firstNodeImageIdx = firstNodeIdx(firstNodeImageIds == imageId);
               secondNodeImageIdx = secondNodeIdx(secondNodeImageIds == imageId);
               firstNodeCoordsInImage = firstNodeCoords(firstNodeImageIds == imageId,:);
               secondNodeCoordsInImage = secondNodeCoords(secondNodeImageIds == imageId,:);
               numberOfFirstNodes = size(firstNodeCoordsInImage,1);

               % Iterate over each seed node and find its adjacentNodes.
               for nodeItr = 1:numberOfFirstNodes
                   %% EDGE INFORMATION EXTRACTION
                   centerArr = [ones(size(secondNodeCoordsInImage,1),1) * firstNodeCoordsInImage(nodeItr,1), ...
                           ones(size(secondNodeCoordsInImage,1),1) * firstNodeCoordsInImage(nodeItr,2)];
                   if strcmp(options.edgeType, 'contour') && ~isempty(mainGraph)
                       % In contour type edges, we work on the first
                       % level to reveal existing neighboring
                       % information. High level compositions are
                       % linked with their adjacency information among
                       % their leaf nodes at the first level.
                       adjacentNodes = zeros(numel(secondNodeImageIdx),1);
                       selfLeafNodes = currentGraphLevel(firstNodeImageIdx(nodeItr)).leafNodes;
                       selfLeafNeighbors = cell2mat(leafNodeAdjArr(selfLeafNodes));
                       % Estimate edge novelty threshold so that we can
                       % connect nodes that provide novel information.
                       secondNodeCount = numel(secondNodeImageIdx);
                       for secNodeItr = 1:secondNodeCount
                          neighborLeaveCount = sum(ismember(currentGraphLevel(secondNodeImageIdx(secNodeItr)).leafNodes, selfLeafNeighbors));
                          if neighborLeaveCount > 0
                             adjacentNodes(secNodeItr) = 1;
                          end
                       end
                       adjacentNodes = find(adjacentNodes);
                   else
                        distances = sqrt(sum((centerArr - secondNodeCoordsInImage).^2, 2));
                        [adjacentNodes] = find(distances <= neighborhood);
                   end

                   % Need to remove the node itself if checking neighbors with same
                   % label id.
                   if node1 == node2
                        adjacentNodes = adjacentNodes(adjacentNodes ~= nodeItr,:);
                   end
                   relativeVector = secondNodeCoordsInImage - centerArr;
                   relativeVector = relativeVector(adjacentNodes,:);
                   samples = [samples; relativeVector];
               end
           end
           
           %% Cluster samples to detect the nodes. 
           % TOCHANGE: This step can be replaced with a different approach.
           classes = assignModes(samples, options);
           
           %% Calculate cluster centers and reassign the samples to clusters.
           % This step is required since we would like to have
           % consistency with our test cases.
           numberOfClusters = max(classes);
           centers = zeros(numberOfClusters,4);
           centers(:,1) = ones(numberOfClusters,1) * node1;
           centers(:,2) = ones(numberOfClusters,1) * node2;
           for centerItr = 1:numberOfClusters
              centers(centerItr,3:4) = mean(samples(classes==centerItr,:),1);
           end
           % Assign centers to modes.
           modes = [modes; centers];

           % Change mode offset and move on.
           modeOffset = modeOffset + max(classes);
           %% In debug mode, write classes to the output as images.
           if options.debug
               distributionImg = zeros(options.maxImageDim*2+1);
               samplesToWrite = floor(samples + options.maxImageDim + 1);
               samplesInd = sub2ind(size(distributionImg), samplesToWrite(:,1), samplesToWrite(:,2));
               distributionImg(samplesInd) = classes;

               % Resize the distribution image so it is of the smallest
               % possible size.
               [posSampleX, posSampleY] = find(distributionImg);
               distributionImg = distributionImg(min(posSampleX):max(posSampleX), min(posSampleY):max(posSampleY));
               if ~exist([options.currentFolder '/debug/' datasetName '/level' num2str(currentLevel) '/pairwise/'], 'dir')
                   mkdir([options.currentFolder '/debug/' datasetName '/level' num2str(currentLevel) '/pairwise/']);
               end
               if ~isempty(distributionImg)
                   imwrite(label2rgb(distributionImg, 'jet', 'k', 'shuffle'), [options.currentFolder '/debug/' datasetName '/level' num2str(currentLevel) '/pairwise/' num2str(node1) '_' num2str(node2) '.png']);
               end
           end
       end
    end
    warning(w);
end
