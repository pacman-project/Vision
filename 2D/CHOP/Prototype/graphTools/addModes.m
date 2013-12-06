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
%> @param options Program options.
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
function [modes, edges] = addModes(nodes, options, currentLevel)
    % Prevent empty cluster warnings in kmeans.
    w = warning('off', 'stats:kmeans:EmptyCluster');
    
    % Calculate edge radius.
    scale = (1/options.scaling)^currentLevel;
    neighborhood = options.edgeRadius * scale;
    
    %% For each node pair, get the 2-D distribution of samples.
    % Right now, we only work with 2-dimensional spatial relations.
    numberOfNodes = size(nodes,1);
    nodeIds = cell2mat(nodes(:,1));
    nodeCoords = cell2mat(nodes(:,2));
    imageIds = cell2mat(nodes(:,3));
    
    modes = [];
    edges = [];
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
               firstNodeCoordsInImage = firstNodeCoords(firstNodeImageIds == imageId,:);
               secondNodeCoordsInImage = secondNodeCoords(secondNodeImageIds == imageId,:);
               numberOfFirstNodes = size(firstNodeCoordsInImage,1);
               
               firstNodesInImage = firstNodeIdx(firstNodeImageIds == imageId,:);
               secondNodesInImage = secondNodeIdx(secondNodeImageIds == imageId,:);
               
               % Iterate over each seed node and find its adjacentNodes.
               for nodeItr = 1:numberOfFirstNodes
                   centerArr = [ones(size(secondNodeCoordsInImage,1),1) * firstNodeCoordsInImage(nodeItr,1), ...
                       ones(size(secondNodeCoordsInImage,1),1) * firstNodeCoordsInImage(nodeItr,2)];
               
                   distances = sqrt(sum((centerArr - secondNodeCoordsInImage).^2, 2));
                   adjacentNodes = find(distances <= neighborhood);
                   adjacentNodes = adjacentNodes(adjacentNodes ~= nodeItr);
                   relativeVector = secondNodeCoordsInImage - centerArr;
                   relativeVector = relativeVector(adjacentNodes,:);
                   samples = [samples; relativeVector];
                   
                   %% If both nodes are the same, we put twice directed edges to make the final graph undirected.
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
           
           %% Cluster samples to detect the nodes. 
           % This step can be replaced with a better, different approach.
           if isempty(samples)
               continue;
           elseif size(samples,1) < options.maximumModes
               % If not enough samples, combine all into single cluster.
               classes = ones(size(samples,1),1);
           else
               % Enough samples, process data.
               classes = mec(samples, 'c', options.maximumModes, 'kmeans_i', 5);
           end
           allPairwiseEdges(:,3) = classes+modeOffset;
           edges = [edges; allPairwiseEdges];
           modeOffset = modeOffset + max(classes);
           
           %% Estimate cluster centers out of classes.
           % These will be used to reconstruct the object in the top-down
           % parsing.
           numberOfClusters = max(classes);
           centers = zeros(numberOfClusters,4);
           centers(:,1) = ones(numberOfClusters,1) * node1;
           centers(:,2) = ones(numberOfClusters,1) * node2;
           for centerItr = 1:numberOfClusters
              centers(centerItr,3:4) = mean(samples(classes==centerItr,:),1);
           end
           % Assign centers to modes.
           modes = [modes; centers];
           
           %% In debug mode, write classes to the output as images.
           if options.debug
               distributionImg = zeros(neighborhood*2+1);
               samplesToWrite = floor(samples + neighborhood + 1);
               samplesInd = sub2ind(size(distributionImg), samplesToWrite(:,1), samplesToWrite(:,2));
               distributionImg(samplesInd) = classes;
               if ~exist([options.currentFolder '/debug'], 'dir')
                   mkdir([options.currentFolder '/debug']);
               end
               if ~exist([options.currentFolder '/debug/level' num2str(currentLevel)], 'dir')
                   mkdir([options.currentFolder '/debug/level' num2str(currentLevel)]);
               end
      %         imwrite(label2rgb(distributionImg, 'jet', 'k', 'shuffle'), [options.currentFolder '/debug/level' currentLevel '/' num2str(node1) '_' num2str(node2) '.png']);
           end
       end
    end
    warning(w);
end

