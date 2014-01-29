%> Name: createEdgesWithLabels
%>
%> Description: Given the node list, the aim of this
%> function is to detect spatial distributions by analyzing 2-D spatial
%> arrangements of node types. For example, if there are 3 distinct
%> configurations of a node type pair, the edges in between are represented 
%> with 3 categories.
%>
%> @param nodes The node list including label ids, positions and image ids
%>      of each node.
%> @param mainGraph The object graphs' data structure.
%> @param leafNodeAdjArr The adjacency array of leaf nodes, in line with
%>      regular edge generation in the hierarchy. Neighboring low-level
%>      realizations in first layer are linked together. In first level,
%> this parameter is empty, and is returned by this function. In upper
%> levels, it is given as a valid (non-empty) parameter.
%> @param options Program options.
%> @param currentLevel The current scene graph level.
%> @param modes If empty, new modes are to be learned. If not, edge labels
%>      will be formed depending on existing modes.
%> 
%> @retval edges Edges are of the form: [ node1, node2, mode, directed;
%>                                        node1, node2, mode, directed;
%>                                      ...]
%> @retval leafNodeAdjArr Valid response in level 1, where leaf nodes are
%> linked if they are actually linked in the object graphs.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 04.12.2013
%> Separate mode learning from this function on 28.01.2014
function [edges, leafNodeAdjArr] = createEdgesWithLabels(nodes, mainGraph, leafNodeAdjArr, options, currentLevel, modes)

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
    edges = zeros(options.maxNumberOfEdges,4);
    edgeOffset = 1;
    
    %% In level 1, we create first level adjacency information between nodes.
    if currentLevel==1 && strcmp(options.edgeType, 'contour')
        leafNodeAdjArr = cell(size(nodes,1),1);
    end
    
    % Read histogram matrix to be used in 'hist' type edge calculations.
    if strcmp(options.property, 'hist')
       load('hMatrix.mat', 'hMatrix'); 
       sizeHMatrix = size(hMatrix,1);
       halfSizeHMatrix = floor(sizeHMatrix/2);
    end
    
    %% If necessary, learn pair-wise node distributions (modes)
    if ~isempty(modes) && numel(modes) >= currentLevel
        currentModes = modes{currentLevel};
    else
        currentModes = [];
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
               %% Verify that the adjacent nodes found by 'contour' info are actually in the neighborhood of seed node.
               distances = sqrt(sum((centerArr - imageCoords).^2, 2));
               adjacentNodes = distances <= neighborhood & adjacentNodes;
               adjacentNodes = find(adjacentNodes);
           else
               % 'distance' type score calculation
               distances = sqrt(sum((centerArr - imageCoords).^2, 2));
               adjacentNodes = find(distances <= neighborhood);
               
               % Looking for closest nodes to connect!
               scores = distances(adjacentNodes);  
           end
           
           %% If no adjacent nodes exist, exit.
           if isempty(adjacentNodes) 
               continue;
           end
           
           %% Get valid adjacent nodes with larger node ids.
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

           %% Create an edge for each seed-adjacent node pair.
           for adjItr = 1:numel(adjacentNodes)
               sample = imageCoords(adjacentNodes(adjItr),:) - imageCoords(nodeItr,:);
               % Decide whether to put a directed link.
               if imageNodeIds(nodeItr) == imageNodeIds(adjacentNodes(adjItr))
                   isDirected = 1;
                   % Check if this edge is on the right side of separating line.
                   if (sample(1) >= 0 && sum(sample) < 0) || (sample(1) < 0 && sum(sample) <= 0)
                       continue;
                   end
                   
                   % If both are at the origin (same spot), this causes
                   % into duplicate links in our final graph. To prevent
                   % this, we only link one with smaller index to the
                   % other one with larger index.
                   if sample(1) == 0 && sample(2) == 0 && nodeItr >= adjacentNodes(adjItr)
                      continue; 
                   end
               else
                   isDirected = 0;
               end

               %% Assign this edge its label based on property to be examined.
               if strcmp(options.property, 'mode')
                    % Estimate mode of this edge.
                    applicableModeIdx = find(currentModes(:,1) == imageNodeIds(nodeItr) & ...
                    currentModes(:,2) == imageNodeIds(adjacentNodes(adjItr)));
                    applicableModes = currentModes(applicableModeIdx,3:4);
                    centerArr2 = [ones(size(applicableModes,1),1) * sample(1), ones(size(applicableModes,1),1) * sample(2)];
                    distances = sqrt(sum((centerArr2 - applicableModes).^2, 2));
                    [~, minDist] = min(distances);
                    
                    % If a valid mode exists, assign its label.
                    if numel(minDist)>0
                        edgeId = applicableModeIdx(minDist(1));
                    else
                        edgeId = 0;
                    end
               elseif strcmp(options.property, 'hist')
                    % Assign the edge a label based on the histogram info
                    normalizedSample = round((double(sample) / neighborhood)*halfSizeHMatrix) + halfSizeHMatrix;
                    normalizedSample(normalizedSample < 1) = 1;
                    normalizedSample(normalizedSample > sizeHMatrix) = sizeHMatrix;
                    edgeId = hMatrix(normalizedSample(1), normalizedSample(2));
               else 
                    % Uniform labeling ('co-occurence' property)
                    edgeId = 1;
               end
               
               %% Create the edge with the estimated label.
               if edgeId ~= 0
                   edges(edgeOffset,:) = [nodeItr + nodeOffset, adjacentNodes(adjItr) + nodeOffset, edgeId, isDirected];
                   edgeOffset = edgeOffset + 1;
               end
           end
        end
        nodeOffset = nodeOffset + numberOfNodes;
    end
    
    % Get nonzero rows of edges
    edges = edges(1:(edgeOffset-1),:);
    
    %% Remove duplicate neighbors from leaf node adjacency list.
    if currentLevel==1 && strcmp(options.edgeType, 'contour')
        leafNodeAdjArr = cellfun(@(x) unique(x), leafNodeAdjArr, 'UniformOutput', false);
    end
end
