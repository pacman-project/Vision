%> Name: createEdgesWithLabels
%>
%> Description: Given the node list, the aim of this
%> function is to detect spatial distributions by analyzing 2-D spatial
%> arrangements of node types. For example, if there are 3 distinct
%> configurations of a node type pair, the edges in between are represented 
%> with 3 categories.
%>
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
function [mainGraph, leafNodeAdjArr] = createEdgesWithLabels(mainGraph, leafNodeAdjArr, options, currentLevelId, modes)
    %% Function initializations, reading data from main graph.
    % Calculate edge radius.
    scale = (1/options.scaling)^(currentLevelId-1);
    neighborhood = fix(options.edgeRadius * scale);
    
    currentLevel = mainGraph{currentLevelId};
    nodeIds = [currentLevel.labelId]';
    nodeCoords = cat(1, currentLevel.position);
    imageIds = [currentLevel.imageId]';
    rfIds = [currentLevel.rfId]';
    isCenter = [currentLevel.isCenter]';
    numberOfNodes = numel(currentLevel);
    
    %% In level 1, we create first level adjacency information between nodes.
    % Else, it is already given.
    if currentLevelId==1 && strcmp(options.edgeType, 'contour')
        leafNodeAdjArr = cell(numberOfNodes,1);
    end
    
    %% Read histogram matrix to be used in 'hist' type edge calculations.
    if strcmp(options.property, 'hist')
       load('hMatrix.mat', 'hMatrix'); 
       sizeHMatrix = size(hMatrix,1);
       halfSizeHMatrix = floor(sizeHMatrix/2);
    end
    
    %% Get relevant pair-wise distributions (nodes)
    if ~isempty(modes) && numel(modes) >= currentLevelId
        currentModes = modes{currentLevelId};
    else
        currentModes = [];
    end
    
    %% Assuming the distributions have already been learned, create edges and assign edge labels accordingly.
    nodeOffset = 0;
    for imageItr = 1:max(imageIds)
        imageNodeIdx = imageIds == imageItr;
        imageRfIds = rfIds(imageNodeIdx,:);

        % If there are no nodes in this image, move on.
        if nnz(imageNodeIdx) == 0
           continue; 
        end

        %% Process each receptive field in this image. 
        % If only receptive field id is zero, that means we do not
        % implement receptive fields, and the graph generation can proceed
        % as it would. If there are other integer ids, basically each of
        % them should be considered separately.
        rfIdSet = unique(imageRfIds);
        for rfId = rfIdSet'
            
            % Prepare receptive-field related data structures.
            rfNodeIdx = find(imageNodeIdx & rfIds == rfId);
            rfIsCenter = find(isCenter(rfNodeIdx));
            rfNodeIds = nodeIds(rfNodeIdx,:);
            rfCoords = nodeCoords(rfNodeIdx,:);
            numberOfRfNodes = size(rfCoords,1);
            
            for nodeItr = rfIsCenter'
               %% Depending on the edge type, find a list of nodes which are neighbors of this node.
               centerArr = repmat(rfCoords(nodeItr,:), size(rfNodeIds,1),1);
               distances = sqrt(sum((centerArr - rfCoords).^2, 2));
               nearbyNodes = distances <= neighborhood;
               
               % Check adjacent nodes to see if they are among nodes
               % introducing enough novelty with respect to each other. 
               % options.edgeNoveltyThr provides a lower bound for the
               % novelty each edge should introduce.
               if currentLevelId > 1
                   rfGraphNodes = currentLevel(rfNodeIdx);
                   centerLeafNodes = currentLevel(rfNodeIdx(nodeItr)).leafNodes;
                   rfLeafNodes = {rfGraphNodes.leafNodes};
                   leafCounts = cellfun(@(x) numel(x), rfLeafNodes);
                   diffLeafCounts = cellfun(@(x) numel(setdiff(x, centerLeafNodes)), rfLeafNodes);
                   novelNodes = (diffLeafCounts >= options.edgeNoveltyThr*leafCounts)';
                   
                   adjacentNodes = novelNodes & nearbyNodes;
                   % Keep intersecting nodes here for contour type edges.
                   if strcmp(options.edgeType, 'contour')
                       
                   end
               else
                   adjacentNodes = nearbyNodes;
               end
               
               if strcmp(options.edgeType, 'contour')
                   % In contour type edges, we work on the first
                   % level to reveal existing neighboring
                   % information. High level compositions are
                   % linked with their adjacency information among
                   % their leaf nodes at the first level.
                   selfLeafNodes = currentLevel(imageNodeIdx(nodeItr)).leafNodes;
                   selfLeafNeighbors = cell2mat(leafNodeAdjArr(selfLeafNodes));
                   scores = zeros(numel(imageNodeIdx),1);

                   %% Score calculation for each edge takes place in here.
                   for secNodeItr = 1:numel(imageNodeIdx)
                      bondStrength = sum(ismember(currentLevel(imageNodeIdx(secNodeItr)).leafNodes, selfLeafNeighbors));
                      maxOverlap = numel(currentLevel(imageNodeIdx(secNodeItr)).leafNodes) * (1-options.edgeNoveltyThr);
                      numberOfOverlapLeaves = sum(ismember(selfLeafNodes, currentLevel(imageNodeIdx(secNodeItr)).leafNodes));
                      if currentLevelId == 1 || ((bondStrength > 0 && numberOfOverlapLeaves <= maxOverlap))
                         scores(secNodeItr) = numberOfOverlapLeaves - bondStrength;
                         adjacentNodes(secNodeItr) = 1;
                      end
                   end
                   %% Verify that the adjacent nodes found by 'contour' info are actually in the neighborhood of seed node.
                   distances = sqrt(sum((centerArr - rfCoords).^2, 2));
                   adjacentNodes = distances <= neighborhood & adjacentNodes;
                   adjacentNodes = find(adjacentNodes);
               else
                   % 'distance' type score calculation
                   distances = sqrt(sum((centerArr - rfCoords).^2, 2));
                   adjacentNodes = distances <= neighborhood & adjacentNodes;
                   adjacentNodes = find(adjacentNodes);

                   % Looking for closest nodes to connect!
                   scores = distances(adjacentNodes);  
               end

               
               % Let's remove center node from its neighbors.
               otherAdjacentNodes = adjacentNodes ~= nodeItr;
               adjacentNodes = adjacentNodes(otherAdjacentNodes);
               scores = scores(otherAdjacentNodes);
               
               % If no adjacent nodes exist, exit.
               if isempty(adjacentNodes)
                   continue;
               end
               
               %% If using receptive fields, there is no limit on the number of neighbors.
               if ~options.useReceptiveField
                   % Get valid adjacent nodes with larger node ids.
                   validAdjacentNodes = rfNodeIds(adjacentNodes) >= rfNodeIds(nodeItr);
                   adjacentNodes = adjacentNodes(validAdjacentNodes);
                   scores = scores(validAdjacentNodes);

                   % Eliminate low-scored adjacency links to keep the graph degree at a constant level.
                   if currentLevelId == 1
                       averageNodeDegree = options.maxNodeDegreeLevel1;
                   else
                       averageNodeDegree = options.maxNodeDegree;
                   end
                   if numel(adjacentNodes)>averageNodeDegree
                        [idx] = getSmallestNElements(scores, averageNodeDegree);
                        adjacentNodes = adjacentNodes(idx);
                   end
               end
               %% For level 1, add adjacency information to the leafNodeAdjArr array
               % Please note that this
               % adjacency information is directed (a->b =/= b->a)
               if strcmp(options.edgeType, 'contour') && currentLevelId == 1
                   leafNodeAdjArr(rfNodeIdx(nodeItr)) = {[leafNodeAdjArr{rfNodeIdx(nodeItr)}; rfNodeIdx(adjacentNodes)]};
               end
               
               %% Create an edge for each center-adjacent node pair.
               edges = zeros(numel(adjacentNodes),4);
               edgeOffset = 1;
               for adjNode = adjacentNodes'
                   sample = rfCoords(adjNode,:) - rfCoords(nodeItr,:);
                   if ~options.useReceptiveField
                       % Decide whether to put a directed link.
                       if rfNodeIds(nodeItr) == rfNodeIds(adjNode)
                           isDirected = 1;
                           % Check if this edge is on the right side of separating line.
                           if (sample(1) >= 0 && sum(sample) < 0) || (sample(1) < 0 && sum(sample) <= 0)
                               continue;
                           end

                           % If both are at the origin (same spot), this causes
                           % into duplicate links in our final graph. To prevent
                           % this, we only link one with smaller index to the
                           % other one with larger index.
                           if sample(1) == 0 && sample(2) == 0 && nodeItr >= adjNode
                              continue; 
                           end
                       else
                           isDirected = 0;
                       end
                   else
                       isDirected = 1;
                   end

                   %% ASSIGN EDGE LABEL BASED ON PROPERTY HERE.
                   if strcmp(options.property, 'mode')
                        % Estimate mode of this edge.
                        applicableModeIdx = find(currentModes(:,1) == rfNodeIds(nodeItr) & ...
                        currentModes(:,2) == rfNodeIds(adjNode));
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
                       edges(edgeOffset,:) = [nodeItr + nodeOffset, adjNode + nodeOffset, edgeId, isDirected];
                       edgeOffset = edgeOffset + 1;
                   end
               end
               %% Assign all edges to the node in the given graph.
               if edgeOffset > 1
                   currentLevel(rfNodeIdx(nodeItr)).adjInfo = edges((1:edgeOffset-1), :);
               end
            end
            nodeOffset = nodeOffset + numberOfRfNodes;
        end
    end
    mainGraph(currentLevelId) = {currentLevel};
    
    %% Remove duplicate neighbors from leaf node adjacency list.
    if currentLevelId==1 && strcmp(options.edgeType, 'contour')
        leafNodeAdjArr = cellfun(@(x) unique(x), leafNodeAdjArr, 'UniformOutput', false);
    end
end
