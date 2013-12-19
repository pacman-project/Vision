%> Name: addModes
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
function [modes, edges] = addModes(nodes, mainGraph, options, currentLevel, datasetName, modes)
    % Prevent empty cluster warnings in kmeans.
    w = warning('off', 'all');
    
    % Calculate edge radius.
    scale = (1/options.scaling)^(currentLevel-1);
    neighborhood = fix(options.edgeRadius * scale);
    
    nodeIds = cell2mat(nodes(:,1));
    nodeCoords = cell2mat(nodes(:,2));
    imageIds = cell2mat(nodes(:,3));
    numberOfNodes = max(nodeIds);
    if strcmp(options.edgeType, 'contour') && ~isempty(mainGraph)
        firstGraphLevel = mainGraph{1};
        currentGraphLevel = mainGraph{currentLevel};
    end
    edges = [];
    
    if isempty(modes)
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
               allPairwiseEdges = [];
               % Find all edges in between.
               for imageId = 1:max(firstNodeImageIds)
                   firstNodeImageIdx = firstNodeIdx(firstNodeImageIds == imageId);
                   secondNodeImageIdx = secondNodeIdx(secondNodeImageIds == imageId);
                   firstNodeCoordsInImage = firstNodeCoords(firstNodeImageIds == imageId,:);
                   secondNodeCoordsInImage = secondNodeCoords(secondNodeImageIds == imageId,:);
                   numberOfFirstNodes = size(firstNodeCoordsInImage,1);

                   firstNodesInImage = firstNodeIdx(firstNodeImageIds == imageId,:);
                   secondNodesInImage = secondNodeIdx(secondNodeImageIds == imageId,:);

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
                           adjacentNodes = [];
                           selfLeafNodes = currentGraphLevel(firstNodeImageIdx(nodeItr)).leafNodes;
                           selfLeafNeighbors = [];
                           for selfLeafItr = 1:numel(selfLeafNodes)
                                tempNeighbors = unique(firstGraphLevel(selfLeafNodes(selfLeafItr)).adjInfo(:,1:2));
                                tempNeighbors = setdiff(tempNeighbors, selfLeafNodes(selfLeafItr));
                                selfLeafNeighbors = [selfLeafNeighbors; tempNeighbors];
                           end
                           for secNodeItr = 1:numel(secondNodeImageIdx)
                              if numel(intersect(selfLeafNeighbors, currentGraphLevel(secondNodeImageIdx(secNodeItr)).leafNodes)) > 0
                                 adjacentNodes = [adjacentNodes; secNodeItr];
                              end
                           end
                       else
                            distances = sqrt(sum((centerArr - secondNodeCoordsInImage).^2, 2));
                            adjacentNodes = find(distances <= neighborhood);    
                       end
                       
                       % Need to remove the node itself if checking neighbors with same
                       % label id.
                       if node1 == node2
                            adjacentNodes = adjacentNodes(adjacentNodes ~= nodeItr);
                       end
                       relativeVector = secondNodeCoordsInImage - centerArr;
                       relativeVector = relativeVector(adjacentNodes,:);
                       samples = [samples; relativeVector];

                       %% If both nodes are the same, we put directed edges.
                       numberOfCurrEdges = size(adjacentNodes,1);

                       % Create edges here.
                       if node1 == node2
                           directedArr = ones(numberOfCurrEdges,1);
                       else
                           directedArr = zeros(numberOfCurrEdges,1);
                       end

                       currEdges = [ones(numberOfCurrEdges,1) * firstNodesInImage(nodeItr), ...
                           secondNodesInImage(adjacentNodes), zeros(numberOfCurrEdges,1), directedArr];
                       if ~isempty(currEdges)
                           allPairwiseEdges = [allPairwiseEdges; currEdges];
                       end
                   end
               end
               % Normalize vectors using max distance
     %          samples = samples / neighborhood;

               % Eliminate some samples so that we only take roughly half of
               % them into account (only when investigating same node type
               % pair).
               if node1 == node2 && ~isempty(samples)
                   validEdgeIdx = samples(:,1) + samples(:,2) >= 0;
                   samples = samples(validEdgeIdx, :);
                   allPairwiseEdges = allPairwiseEdges(validEdgeIdx, :);
               end

               %% Cluster samples to detect the nodes. 
               % This step can be replaced with a better, different approach.
               if isempty(samples)
                   continue;
               elseif size(samples,1) < options.maximumModes
                   if size(samples,1) > 1
                        classes = mec(samples, 'c', size(samples,1), 'kmeans_i', 5);
                   else
                        classes = 1;
                   end
               else
                   % Enough samples, process data.
                   classes = mec(samples, 'c', options.maximumModes, 'kmeans_i', 5);
               end
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
               
               % Reassign the samples to the clusters.
               for sampleItr = 1:size(samples,1)
                  sampleArr = [ones(size(centers,1),1)* samples(sampleItr,1), ones(size(centers,1),1)* samples(sampleItr,2)];
                  distances = sqrt(sum((sampleArr - centers(:,3:4)).^2, 2));
                  [~, minDist] = min(distances);
                  if numel(minDist)>0
                     classes(sampleItr) = minDist(1);
                  end
               end
               
               % Assign edges, change mode offset and move on.
               allPairwiseEdges(:,3) = classes+modeOffset;
               edges = [edges; allPairwiseEdges];
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
                   imwrite(label2rgb(distributionImg, 'jet', 'k', 'shuffle'), [options.currentFolder '/debug/' datasetName '/level' num2str(currentLevel) '/pairwise/' num2str(node1) '_' num2str(node2) '.png']);
               end
           end
        end
        warning(w);

    else
        %% Assuming the distributions have already been learned, create edges and assign edge labels accordingly.
        currentModes = modes{currentLevel};
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
                   adjacentNodes = [];
                   selfLeafNodes = currentGraphLevel(imageNodeIdx(nodeItr)).leafNodes;
                   selfLeafNeighbors = [];
                   for selfLeafItr = 1:numel(selfLeafNodes)
                        tempNeighbors = unique(firstGraphLevel(selfLeafNodes(selfLeafItr)).adjInfo(:,1:2));
                        tempNeighbors = setdiff(tempNeighbors, selfLeafNodes(selfLeafItr));
                        selfLeafNeighbors = [selfLeafNeighbors; tempNeighbors];
                   end
                   for secNodeItr = 1:numel(imageNodeIdx)
                      if numel(intersect(selfLeafNeighbors, currentGraphLevel(imageNodeIdx(secNodeItr)).leafNodes)) > 0
                         adjacentNodes = [adjacentNodes; secNodeItr];
                      end
                   end
               else
                    distances = sqrt(sum((centerArr - imageCoords).^2, 2));
                    adjacentNodes = find(distances <= neighborhood);  
               end
               adjacentNodes = adjacentNodes(adjacentNodes ~= nodeItr);
               adjacentNodes = adjacentNodes(imageNodeIds(adjacentNodes) >= imageNodeIds(nodeItr));
               
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
    end
end

