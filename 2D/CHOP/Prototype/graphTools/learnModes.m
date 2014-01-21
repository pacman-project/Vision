%> Name: learnModes
%>
%> Description: Given the node list of all images, this function learns the
%> modes by clustering pairwise relative positions in 2D space. 
%> @param nodes The node list including label ids, positions and image ids
%> of each node.
%> @param mainGraph The object graphs' data structure.
%> @param leafNodeAdjArr Leaf node adjacency array, which is basically an
%> adjacency list of every leaf node that is connected because they are in
%> the neighborhood of each other.
%> @param options Program options.
%> @param currentLevel The current scene graph level.
%> @param datasetName Name of the dataset.
%> 
%> @retval modes The mode list representing edge categories.
%> @retval edges Edges are of the form: [ node1, node2, mode, directed;
%>                                        node1, node2, mode, directed;
%>                                      ...]
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 21.01.2014
function [modes] = learnModes(nodes, mainGraph, leafNodeAdjArr, options, currentLevel, datasetName)
    % Prevent empty cluster warnings in kmeans.
    w = warning('off', 'all');
    
    % Calculate edge radius.
    scale = (1/options.scaling)^(currentLevel-1);
    neighborhood = fix(options.edgeRadius * scale);
    
    % Set initial data structures for processing 
    modes = [];
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
                       
                       %% Verify that the adjacent nodes found by 'contour' info are actually in the neighborhood of seed node.
                       distances = sqrt(sum((centerArr - imageCoords).^2, 2));
                       adjacentNodes = distances <= neighborhood & adjacentNodes;
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
                   % Get rid of half where 
                   % relativeVector(1) + relativeVector(2) < 0, if node1 == node2.
                   if node1 == node2
                        relativeVector = relativeVector(sum(relativeVector,2) >= 0,:);
                   end
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
               distributionImgDim = max(max(samples)) * 2;
               
               % If no samples are to be written, move on.
               if numel(samplesToWrite) < 1
                   continue;
               end
                   
               samplesInd = sub2ind(size(distributionImg), samplesToWrite(:,1), samplesToWrite(:,2));
               distributionImg(samplesInd) = classes;

               % Resize the distribution image so it is of the smallest
               % possible size.
               [posSampleX, posSampleY] = find(distributionImg);
               minX = min(posSampleX);
               maxX = max(abs(minX), max(posSampleX));
               minX = maxX - distributionImgDim;
               minY = min(posSampleY);
               maxY = max(abs(minY), max(posSampleY));
               minY = maxY - distributionImgDim;
               distributionImg = distributionImg(minX:maxX, minY:maxY);
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