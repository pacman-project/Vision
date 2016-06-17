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
function [mainGraph] = createEdgesWithLabels(mainGraph, options, currentLevelId, edgeNoveltyThr)
    %% Function initializations, reading data from main graph.
    % Calculate edge radius.
    currentLevel = mainGraph{currentLevelId};
    firstLevel = mainGraph{1};
    firstLevelAdjInfo = {firstLevel.adjInfo};
    nodeIds = [currentLevel.labelId]';
    nodeCoords = cat(1, currentLevel.position);
    imageIds = [currentLevel.imageId]';
    edgeIdMatrix = options.edgeIdMatrix;
    halfMatrixSize = (options.receptiveFieldSize+1)/2;
    incHalfMatrixSize = halfMatrixSize - 1;
    if ismember(currentLevelId, options.smallRFLayers)
        rfRadius = ceil(halfMatrixSize/2);
    else
        rfRadius = halfMatrixSize;
    end
    circularRF = options.circularRF;
    matrixSize = [options.receptiveFieldSize, options.receptiveFieldSize];
    receptiveFieldSize = options.receptiveFieldSize;
    edgeType = options.edgeType;
    minimalEdgeCount = options.minimalEdgeCount;
    imagesPerSet = 10;
    maxFirstLevelNodeDegree = options.maxFirstLevelNodeDegree;
    
    %% Create elliptical RFs for layer 1.
    if currentLevelId == 1
       validMatrices = cell(options.numberOfGaborFilters, 1); 
       smallSize = 1;
       mainSize = 4;
       dummyMask = zeros(receptiveFieldSize) > 0;
       dummyMask(halfMatrixSize-mainSize:halfMatrixSize+mainSize, halfMatrixSize-smallSize:halfMatrixSize+smallSize) = 1;
       angleStep = 180 / options.numberOfGaborFilters;
       for filterItr = 1:options.numberOfGaborFilters;
          rotatedMask = imrotate(dummyMask, -(filterItr-1) * angleStep, 'nearest', 'crop');
          validMatrices{filterItr} = rotatedMask;
       end
        secValidMatrices = validMatrices;
       % Create secondary (thinner) validation matrices.
%        secValidMatrices = cell(options.numberOfGaborFilters, 1);
%        dummyMask = zeros(receptiveFieldSize) > 0;
%        dummyMask(halfMatrixSize-mainSize:halfMatrixSize+mainSize, halfMatrixSize) = 1;
%        for filterItr = 1:options.numberOfGaborFilters;
%           rotatedMask = imrotate(dummyMask, -(filterItr-1) * angleStep, 'nearest', 'crop');
%           secValidMatrices{filterItr} = rotatedMask;
%        end
    else
        validMatrices = [];
        secValidMatrices = [];
    end
    
    %% Program options into variables.
    edgeNoveltyThr = 1-edgeNoveltyThr;
    maxNodeDegree = options.maxNodeDegree;
    
    %% Put each image's node set into a different bin.
    numberOfImages = double(max(imageIds));
    imageNodeIdxSets = cell(numberOfImages, 1);
    imageNodeOffsets = zeros(numberOfImages, 1, 'int32');
    imageGraphNodeSets = cell(numberOfImages, 1);
    imageNodeIdArr = cell(numberOfImages,1);
    imageNodeCoordArr = cell(numberOfImages,1);
    nodeOffset = 0;
    for imageItr = 1:max(imageIds)
       imageNodeIdx = find(imageIds == imageItr)';
       imageNodeIdxSets(imageItr) = {imageNodeIdx};
       imageGraphNodeSets(imageItr) = {currentLevel(imageIds == imageItr)};
       imageNodeIdArr(imageItr) = {nodeIds(imageNodeIdx)};
       imageNodeCoordArr(imageItr) = {nodeCoords(imageNodeIdx,:)};
       imageNodeOffsets(imageItr) = nodeOffset;
       nodeOffset = nodeOffset + nnz(imageNodeIdx);
    end
    
    % Put images in sets, for better parallelization
    if numberOfImages > imagesPerSet
         finalSet = rem(numberOfImages, imagesPerSet);
         numberOfSets = floor(numberOfImages/imagesPerSet);
         sets = repmat(imagesPerSet, numberOfSets, 1);
         if finalSet > 0
              sets = [sets; finalSet];
         end
         sets = mat2cell(1:numberOfImages, 1, sets);
    else
         sets = {1:numberOfImages};
         numberOfSets = 1;
    end
    
    % Create data structures for parallel processing.
    setNodeIdxSets = cell(numberOfSets,1);
    setGraphNodeSets = cell(numberOfSets,1);
    setNodeIdArr = cell(numberOfSets,1);
    setNodeCoordArr = cell(numberOfSets,1);
    for setItr = 1:numberOfSets
         imageIdx = sets{setItr};
         setNodeIdxSets{setItr} = imageNodeIdxSets(imageIdx);
         setGraphNodeSets{setItr} = imageGraphNodeSets(imageIdx);
         setNodeIdArr{setItr} = imageNodeIdArr(imageIdx);
         setNodeCoordArr{setItr} = imageNodeCoordArr(imageIdx);
    end
    
    emptyArr = zeros(0, 'int32');
    dummyOneArrs = cell(maxNodeDegree * 2, 1);
    for itr = 1:(maxNodeDegree*2)
       dummyOneArrs(itr) = {ones(itr,1, 'int32')};
    end
    
    %% Process each set separately (and in parallel)
    for setItr = 1:numberOfSets
         imageIdx = sets{setItr};
         imageNodeIdxSets = setNodeIdxSets{setItr};
         imageGraphNodeSets = setGraphNodeSets{setItr};
         imageNodeIdArr = setNodeIdArr{setItr} ;
         imageNodeCoordArr = setNodeCoordArr{setItr};
 %        disp(['Processing set ' num2str(setItr) ', which has ' num2str(numel(imageIdx)) ' images.']);
         
         % For every image in this set, find edges and save them.
         for imageItr = 1:numel(imageIdx)
             imageNodeOffset = imageNodeOffsets(imageIdx(imageItr));
             imageNodeIdx = imageNodeIdxSets{imageItr};
             numberOfNodes = numel(imageNodeIdx);

             % If there are no nodes in this image, move on.
             if numberOfNodes == 0
                continue; 
             end

             % Get data structures containing information about the nodes in this image.
             curNodeIds = imageNodeIdArr{imageItr};
             curNodeCoords = imageNodeCoordArr{imageItr};
             curGraphNodes = imageGraphNodeSets{imageItr};
             curLeafNodes = {curGraphNodes.leafNodes};
             maxSharedLeafNodes = cellfun(@(x) numel(x), curLeafNodes) * edgeNoveltyThr;
             curAdjacentNodes = cell(numberOfNodes,1);

             % If circular RF is desired, we pre-calculate distances.
             if circularRF && numberOfNodes > 1
                  allDistances = squareform(pdist(single(curNodeCoords)));
             else
                  allDistances = 0;
             end
             
             %% Find all edges within this image.
             numberAdjArr = zeros(numberOfNodes,1);
             for nodeItr = 1:numberOfNodes
                centerArr = repmat(curNodeCoords(nodeItr,:), numberOfNodes,1);
                distances = allDistances(:, nodeItr);
                
                %% If the RF is circular, we need to eliminate adjacent nodes further.
                if circularRF
                     adjacentNodes = distances < rfRadius;
                else
                     adjacentNodes = curNodeCoords(:,1) > (centerArr(:,1) - rfRadius) & ...
                          curNodeCoords(:,1) < (centerArr(:,1) + rfRadius) & ...
                          curNodeCoords(:,2) > (centerArr(:,2) - rfRadius) & ...
                          curNodeCoords(:,2) < (centerArr(:,2) + rfRadius);
                end
                adjacentNodes(nodeItr) = 0;
                adjacentNodes = find(adjacentNodes);
                
                %% A further check is done to ensure boundary/surface continuity.
                if strcmp(edgeType, 'continuity') && currentLevelId > 1
                    centerLeafEdges = cat(1, firstLevelAdjInfo{curLeafNodes{nodeItr}}); %#ok<PFBNS>
                    if isempty(centerLeafEdges)
                         adjacentNodes = [];
                    else
                         centerAdjNodes = sort(centerLeafEdges(:,2));
                         adjLeafNodes = curLeafNodes(adjacentNodes);
                         validAdjacentNodes = cellfun(@(x) nnz(ismembc(x, centerAdjNodes)), adjLeafNodes) > 0;
                         if nnz(validAdjacentNodes > 1)
                            adjacentNodes = adjacentNodes(validAdjacentNodes);
                         end
                    end
                end
                
                %% Check for edge novelty.
                if currentLevelId > 1 
                    if edgeNoveltyThr < 1
                        commonLeafCounts = cellfun(@(x) sum(ismembc(x, curLeafNodes{nodeItr})), curLeafNodes(adjacentNodes));
                        novelNodes = (commonLeafCounts <= maxSharedLeafNodes(adjacentNodes))';
                        if nnz(novelNodes) > 1
                            adjacentNodes = adjacentNodes(novelNodes);
                        end
                    end
                else
                    % For layer 1, we check for contour connectivity when
                    % forming an edge between two nodes.
                    centerRF = validMatrices{curNodeIds(nodeItr)};
                    adjacentNodeCoords = (curNodeCoords(adjacentNodes,:) - centerArr(adjacentNodes,:)) + halfMatrixSize;
                    centerCoords = halfMatrixSize - (adjacentNodeCoords - halfMatrixSize); 
                    validAdjacentNodes = ones(numel(adjacentNodes),1) > 0;
                    shareRFCounts = zeros(numel(adjacentNodes),1);
                    for adjItr = 1:numel(adjacentNodes)
                        tempRF = centerRF(max(1, adjacentNodeCoords(adjItr,1) - incHalfMatrixSize):...
                            min(receptiveFieldSize, adjacentNodeCoords(adjItr,1) + incHalfMatrixSize), ...
                            max(1, adjacentNodeCoords(adjItr,2) - incHalfMatrixSize): ...
                            min(receptiveFieldSize, adjacentNodeCoords(adjItr,2) + incHalfMatrixSize));
                        
                        % Now, we reverse the situation and get the part of
                        % the RF that is adjacent-node based.
                        adjRF = secValidMatrices{curNodeIds(adjacentNodes(adjItr))};
                        tempRF = tempRF & adjRF(max(1, centerCoords(adjItr,1) - incHalfMatrixSize):...
                            min(receptiveFieldSize, centerCoords(adjItr,1) + incHalfMatrixSize), ...
                            max(1, centerCoords(adjItr,2) - incHalfMatrixSize): ...
                            min(receptiveFieldSize, centerCoords(adjItr,2) + incHalfMatrixSize));
                        
                        % If both rfs do not match, we move on.
                        if nnz(tempRF) < 2
                            validAdjacentNodes(adjItr) = 0;
                        else
                            shareRFCounts(adjItr) = nnz(tempRF);
                        end
                    end
                    
                    % If this node has too many adjacent node pairs, we
                    % simply don't consider unlikely ones.
                    if nnz(shareRFCounts) > maxFirstLevelNodeDegree
                        [~, idx] = sort(shareRFCounts, 'descend');
                        adjacentNodes = adjacentNodes(sort(idx(1:maxFirstLevelNodeDegree)));
                    else
                        adjacentNodes = adjacentNodes(validAdjacentNodes);
                    end
                end

                %% Now, we'll find a canonical edge representation for this node. 
                % We'll pick adjacent nodes one by one, based on their
                % contribution (addition of novel nodes), and their
                % distance. Ideally, we want to pick closer nodes with high
                % contributions.
                if currentLevelId > 1
                     curNodeList = curLeafNodes{nodeItr};
                     
                     
                     if numel(adjacentNodes) > maxNodeDegree || minimalEdgeCount
                         % Alternative 1.
                         novelNodeCounts = cellfun(@(x) nnz(~ismembc(x, curNodeList)), curLeafNodes(adjacentNodes))';
                         if nnz(novelNodeCounts) > 0
                             sortArr = [-novelNodeCounts, distances(adjacentNodes), single(curNodeIds(adjacentNodes))];
                             [~, idx] = sortrows(sortArr);
                             adjacentNodes = adjacentNodes(sort(idx(1:maxNodeDegree)));
                         else
                             adjacentNodes = [];
                         end
                         
                         % Alternative 2.
%                            selectedNodeCount = 0;
%                           pickedAdjacentNodes = [];
%                           while selectedNodeCount < maxNodeDegree && ~isempty(adjacentNodes)
%                               % Count the number of novel nodes for every adjacent
%                               % node (which wasn't picked).
%                               novelNodeCounts = cellfun(@(x) nnz(~ismembc(x, curNodeList)), curLeafNodes(adjacentNodes))';
%                               if nnz(novelNodeCounts) == 0
%                                    break;
%                               end
%                               
%                               % Eliminate instances that have already been
%                               % used for coverage.
%                               validArr = novelNodeCounts > 0;
% 
%                               % Now, we pick next best adjacent node.
%                       %        sortArr = [-novelNodeCounts, distances(adjacentNodes)];
%                               sortArr = [-novelNodeCounts, curNodeIds(adjacentNodes)];
%                               [~, idx] = sortrows(sortArr);
%                               pickedAdjacentNodes = cat(1, pickedAdjacentNodes, adjacentNodes(idx(1)));
%                               validArr(idx(1)) = 0;
%                               curNodeList = fastsortedunique(sort(cat(2, curNodeList, curLeafNodes{adjacentNodes(idx(1))})));
%                               selectedNodeCount = selectedNodeCount + 1;
%                               
%                               % Update lists.
%                               adjacentNodes = adjacentNodes(validArr);
%                           end
%                           adjacentNodes = pickedAdjacentNodes;
                     end
                else
                   %% Eliminate far away nodes if this one has too many neighbors
                   % Calculate scores (distances).
                   scores = distances(adjacentNodes);

                   % Eliminate nodes having lower scores.
                  if numel(adjacentNodes)>maxNodeDegree
                        [idx] = getLargestNElements(scores, maxNodeDegree);
                        adjacentNodes = adjacentNodes(idx);
                  end
                end

                %% Assign final adjacent nodes.
                numberOfAdjacentNodes = numel(adjacentNodes);
                numberAdjArr(nodeItr) = numberOfAdjacentNodes;
                if numberOfAdjacentNodes > 0
                    curAdjacentNodes(nodeItr) = {adjacentNodes}; 
                else
                    curAdjacentNodes(nodeItr) = {emptyArr};
                end
             end
             
             % Trial
             initialsArr = zeros(sum(numberAdjArr), 1, 'int32');
             offset = 1;
             for itr = 1:numberOfNodes
                 initialsArr(offset:(offset+numberAdjArr(itr)-1)) = itr;
                 offset = offset + numberAdjArr(itr);
             end
             
             % Get rid of empty entries in curAdjacentNodes.
             nonemptyCurAdjacentNodeIdx = numberAdjArr > 0;
             curAdjacentNodes = curAdjacentNodes(nonemptyCurAdjacentNodeIdx);
             

             % Obtain edges and count them.
             allEdges = cat(1, curAdjacentNodes{:});
             
             % Add initial nodes as well.
             allEdges = cat(2, initialsArr, allEdges);
             
             numberOfAllEdges = size(allEdges,1);

             if numberOfAllEdges == 0
                continue;
             end

             % Collect info to form the edges.
             node1Labels = curNodeIds(allEdges(:,1));
             node2Labels = curNodeIds(allEdges(:,2));
             node1Coords = curNodeCoords(allEdges(:,1),:);
             node2Coords = curNodeCoords(allEdges(:,2),:);
             edgeCoords = node2Coords - node1Coords;

             % Remove edges between overlapping (same position) nodes having same id.
             labelEqualityArr = node1Labels == node2Labels;
             validEdges = ~(labelEqualityArr & (edgeCoords(:,1) == 0 & edgeCoords(:,2) == 0));

             % If receptive fields are used, every edge is directed.
             directedArr = ones(nnz(validEdges),1, 'int32');

             % Update data structures based on removed edges.
             allEdges = allEdges(validEdges,:);
             edgeCoords = double(edgeCoords(validEdges,:));
             normalizedEdgeCoords = edgeCoords + halfMatrixSize;

             % Double check in order not to go out of bounds.
             normalizedEdgeCoords(normalizedEdgeCoords < 1) = 1;
             normalizedEdgeCoords(normalizedEdgeCoords > matrixSize(1)) = matrixSize(1);

             %Find edge labels.
             matrixInd = sub2ind(matrixSize, normalizedEdgeCoords(:,1), normalizedEdgeCoords(:,2));
             edgeIds = edgeIdMatrix(matrixInd);
             edges = [allEdges(:,1:2) + imageNodeOffset, edgeIds, directedArr];

             % Due to some approximations in neighborhood calculations, some edgeIds might 
             % have been assigned as 0. Eliminate such cases.
             edges = edges(edgeIds>0,:);

            %% Assign all edges to their respective nodes in the final graph.
            if ~isempty(edges)
                tempEdgeArr = cell(numberOfNodes,1);
                for nodeItr = 1:numberOfNodes
                    edgeIdx = edges(:,1) == (nodeItr + imageNodeOffset);
                    if nnz(edgeIdx) > 0
                        tempEdgeArr{nodeItr} = edges(edgeIdx,:);
                    else
                        tempEdgeArr{nodeItr} = emptyArr;
                    end
                end
                [curGraphNodes.adjInfo] = deal(tempEdgeArr{:});
            end
            imageGraphNodeSets(imageItr) = {curGraphNodes};
         end
         setGraphNodeSets{setItr} = imageGraphNodeSets;
    end
    currentLevel = cat(1, setGraphNodeSets{:});
    currentLevel = [currentLevel{:}];
    mainGraph(currentLevelId) = {currentLevel};
    clearvars -except mainGraph
end
