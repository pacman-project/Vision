%> Name: getEdges
%>
%> Description: Given the nodes with their indices and center coordinates,
%> this function forms the edge structure depending on the geometric property
%> to be examined.
%>
%> @param nodes The nodes extracted from the image in the following format:
%>      [ index x y;
%>        index x y;
%>          ...
%>        index x y]
%> @param options Program options.
%> @param currentLevel The level for which graph is extracted. Needed since
%> the local neighborhood is affected by level.
%>
%> @retval edges Edges belonging to nodes.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 20.11.2013
%> Comment changes on 20.01.2014
function [edges, leafNodeAdjArr] = getEdges(nodes, mainGraph, leafNodeAdjArr, options, currentLevel, datasetName)

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
               % We want to connect distant nodes here.
               scores = -distances(adjacentNodes);  
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

%% Helper function to fetch smallest n elements from vector vect.
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