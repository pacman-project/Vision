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
function [currentLevel] = createEdgesWithLabels(currentLevel, firstLevelAdjNodes, level1Coords, options, currentLevelId, edgeNoveltyThr)
    %% Function initializations, reading data from main graph.
    % Calculate edge radius.
    nodeIds = [currentLevel.labelId]';
    nodeCoords = cat(1, currentLevel.position);
    nodePreciseCoords = cat(1, currentLevel.precisePosition);
    imageIds = [currentLevel.imageId]';
    edgeIdMatrix = options.edgeIdMatrix;
    halfMatrixSize = (options.receptiveFieldSize+1)/2;
    smallHalfMatrixSize = (options.smallReceptiveFieldSize+1)/2;
    if ismember(currentLevelId, options.smallRFLayers)
        rfRadius = smallHalfMatrixSize;
    else
        rfRadius = halfMatrixSize;
    end
    
    if currentLevelId > 1
        firstLevelAdjNodes0 = firstLevelAdjNodes{1};
        firstLevelAdjNodes1 = firstLevelAdjNodes{2};
        firstLevelAdjNodes2 = firstLevelAdjNodes{3};
    else
        firstLevelAdjNodes0 = [];
        firstLevelAdjNodes1 = [];
        firstLevelAdjNodes2 = [];
    end
    
    %% Find pooled positions for layer 1 nodes.
    if strcmp(options.filterType, 'gabor')
         stride = options.gabor.stride;
    else
         stride = options.auto.stride;
    end
    poolDim = options.poolDim;
    
    % Calculate pool factor.
    if currentLevelId> 1
        poolFactor = nnz(~ismembc(2:currentLevelId, options.noPoolingLayers));
    else
        poolFactor = 0;
    end
    
    % Allocate space for counts for every center, and find number of
    % accessible leaf nodes for each receptive field.
    if ~isempty(level1Coords)
         level1CoordsPooled = calculatePooledPositions(level1Coords(:,2:3), poolFactor, poolDim, stride);
    else
         level1CoordsPooled = [];
    end
    
    %% Determine radius and other params.
    % If no pooling has been performed at this layer, and previous layer
    % has small RF, we have a minimum RF limit. 
    if ismember(currentLevelId, options.noPoolingLayers) && ismember(currentLevelId - 1, options.smallRFLayers) && currentLevelId < 7
%        minRFRadius = ceil(halfMatrixSize/2);
        minRFRadius = 1;
    else
        minRFRadius = max(1, min(3, 4-currentLevelId));
%        minRFRadius = 1;
    end
    
    circularRF = options.circularRF;
    matrixSize = [options.receptiveFieldSize, options.receptiveFieldSize];
    edgeType = options.edgeType;
    minimalEdgeCount = options.minimalEdgeCount;
    imagesPerSet = 10;
    
    %% Create elliptical RFs for layer 1.
    if currentLevelId == 1
       validMatrices = cell(options.numberOfGaborFilters, 1); 
  %     smallSize = 1;
 %      mainSize = 4;
       dummyMask = options.filters{1} > 0.01;
       dummyMask = imdilate(dummyMask, [0 1 0; 0 1 0; 0 1 0] > 0, 'full');
       dummyMask = imdilate(dummyMask, [0 1 0; 0 1 0; 0 1 0] > 0, 'full');
       dummyMask = imdilate(dummyMask, [0 1 0; 0 1 0; 0 1 0] > 0, 'full');
       dummyMask = [zeros(size(dummyMask,1),1) > 0, dummyMask, zeros(size(dummyMask,1),1) > 0];
       dummyMask = [zeros(size(dummyMask,1),1) > 0, dummyMask, zeros(size(dummyMask,1),1) > 0];
       dummyMask = [zeros(size(dummyMask,1),1) > 0, dummyMask, zeros(size(dummyMask,1),1) > 0];
       validMatrixSize = size(dummyMask,1);
       validMatrixHalfSize = ceil(size(dummyMask,1) / 2);
       incHalfMatrixSize = validMatrixHalfSize - 1;
  %     dummyMask = zeros(receptiveFieldSize) > 0;
  %     dummyMask(halfMatrixSize-mainSize:halfMatrixSize+mainSize, halfMatrixSize-smallSize:halfMatrixSize+smallSize) = 1;
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
        validMatrixHalfSize = [];
        incHalfMatrixSize = [];
        validMatrixSize = [];
    end
    
    %% Program options into variables.
    edgeShareabilityThr = 1 - edgeNoveltyThr;
    maxShareabilityThr = 1 - options.minEdgeNoveltyThr;
    if currentLevelId == 1
        maxNodeDegree = options.maxFirstLevelNodeDegree;
    else
        maxNodeDegree = options.maxNodeDegree;
    end
    
    %% Put each image's node set into a different bin.
    numberOfImages = double(max(imageIds));
    imageNodeIdxSets = cell(numberOfImages, 1);
    imageNodeOffsets = zeros(numberOfImages, 1, 'int32');
    imageGraphNodeSets = cell(numberOfImages, 1);
    imageNodeIdArr = cell(numberOfImages,1);
    imageNodeCoordArr = cell(numberOfImages,1);
    imageNodePreciseCoordArr = cell(numberOfImages,1);
    nodeOffset = 0;
    for imageItr = 1:max(imageIds)
       imageNodeIdx = find(imageIds == imageItr)';
       imageNodeIdxSets(imageItr) = {imageNodeIdx};
       imageGraphNodeSets(imageItr) = {currentLevel(imageIds == imageItr)};
       imageNodeIdArr(imageItr) = {nodeIds(imageNodeIdx)};
       imageNodeCoordArr(imageItr) = {nodeCoords(imageNodeIdx,:)};
       imageNodePreciseCoordArr(imageItr) = {nodePreciseCoords(imageNodeIdx,:)};
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
              numberOfSets = numberOfSets + 1;
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
    setNodePreciseCoordArr = cell(numberOfSets,1);
    setEdgeDiagnostics = cell(numberOfSets,1);
    for setItr = 1:numberOfSets
         imageIdx = sets{setItr};
         setNodeIdxSets{setItr} = imageNodeIdxSets(imageIdx);
         setGraphNodeSets{setItr} = imageGraphNodeSets(imageIdx);
         setNodeIdArr{setItr} = imageNodeIdArr(imageIdx);
         setNodeCoordArr{setItr} = imageNodeCoordArr(imageIdx);
         setNodePreciseCoordArr{setItr} = imageNodePreciseCoordArr(imageIdx);
    end
    
    emptyArr = zeros(0, 'int32');
    dummyOneArrs = cell(maxNodeDegree * 2, 1);
    for itr = 1:(maxNodeDegree*2)
       dummyOneArrs(itr) = {ones(itr,1, 'int32')};
    end
    
    %% Process each set separately (and in parallel)
    parfor setItr = 1:numberOfSets
         imageIdx = sets{setItr};
         imageNodeIdxSets = setNodeIdxSets{setItr};
         imageGraphNodeSets = setGraphNodeSets{setItr};
         imageNodeIdArr = setNodeIdArr{setItr} ;
         imageNodeCoordArr = setNodeCoordArr{setItr};
         imageNodePreciseCoordArr = setNodePreciseCoordArr{setItr};
         edgeDiagnostics = cell(numel(imageIdx),1);
 %        disp(['Processing set ' num2str(setItr) ', which has ' num2str(numel(imageIdx)) ' images.']);
         removedEdgeCounts = 0;
         
         % For every image in this set, find edges and save them.
         for imageItr = 1:numel(imageIdx)
             imageNodeOffset = imageNodeOffsets(imageIdx(imageItr));
             imageNodeIdx = imageNodeIdxSets{imageItr};
             numberOfNodes = numel(imageNodeIdx);

             % If there are no nodes in this image, move on.
             if numberOfNodes == 0
                continue; 
             end
             
             % Allocate space for edge diagnostics.
             imageEdgeDiagnostics = zeros(numberOfNodes, 1, 'uint8');

             % Get data structures containing information about the nodes in this image.
             curNodeIds = imageNodeIdArr{imageItr};
             curNodeCoords = imageNodeCoordArr{imageItr};
             curNodePreciseCoords = imageNodePreciseCoordArr{imageItr};
             curGraphNodes = imageGraphNodeSets{imageItr};
             curLeafNodes = {curGraphNodes.leafNodes};
             curChildren = {curGraphNodes.children};
             if currentLevelId > 1
                 curChildrenSorted = cellfun(@(x) sort(x), curChildren, 'UniformOutput', false);
                 curCenterNodes = cellfun(@(x) x(1), curChildren);
             end
             allowedSharedLeafNodes = cellfun(@(x) numel(x), curLeafNodes) * edgeShareabilityThr;
             maxSharedLeafNodes = cellfun(@(x) numel(x), curLeafNodes) * maxShareabilityThr;
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
                centerPreciseArr = repmat(curNodePreciseCoords(nodeItr,:), numberOfNodes,1);
                distances = allDistances(:, nodeItr);
                
                % Save info, if there are no other nodes in our RF.
                if numberOfNodes == 1
                    imageEdgeDiagnostics(nodeItr) = 1;
                    continue;
                end
                
                %% 1. First, we find all the adjacent nodes in our receptive field.
                if circularRF
                     adjacentNodes = distances < rfRadius & distances >= minRFRadius;
                else
                     adjacentNodes = curNodeCoords(:,1) > (centerArr(:,1) - rfRadius) & ...
                          curNodeCoords(:,1) <= (centerArr(:,1) - minRFRadius) & ...
                          curNodeCoords(:,1) < (centerArr(:,1) + rfRadius) & ...
                          curNodeCoords(:,1) >= (centerArr(:,1) + minRFRadius) & ...
                          curNodeCoords(:,2) > (centerArr(:,2) - rfRadius) & ...
                          curNodeCoords(:,2) <= (centerArr(:,2) - minRFRadius) & ...
                          curNodeCoords(:,2) < (centerArr(:,2) + rfRadius) & ...
                          curNodeCoords(:,2) >= (centerArr(:,2) + minRFRadius);
                end
                adjacentNodes(nodeItr) = 0;
                adjacentNodes = find(adjacentNodes);
                
                % Save info, if there are no other nodes in our RF.
                if numel(adjacentNodes) == 0
                    imageEdgeDiagnostics(nodeItr) = 2;
                    continue;
                end
                 
                %% A further check is done to ensure boundary/surface continuity.
                adjacentNodeNeigh = ones(numel(adjacentNodes),1,'uint8') * 255;
                if strcmp(edgeType, 'continuity') 
                    if currentLevelId > 1
                        for neighItr = 1:max(1, min(3, currentLevelId-1))
                            % Obtain the correct set of adjacency
                            % information first.
                            if neighItr == 1
                                centerAdjNodes = firstLevelAdjNodes0(curLeafNodes{nodeItr}, :);
                            elseif neighItr == 2
                                centerAdjNodes = firstLevelAdjNodes1(curLeafNodes{nodeItr}, :);
                            else
                                centerAdjNodes = firstLevelAdjNodes2(curLeafNodes{nodeItr}, :);
                            end
                            centerAdjNodes = centerAdjNodes(:)';
                            centerAdjNodes = centerAdjNodes(centerAdjNodes>0);
                            
                            % Check if center's adjacent leaf nodes exist in any
                            % of our neighbors.
                            if isempty(centerAdjNodes)
                                 break;
                            else
                                 centerAdjNodes = fastsortedunique(sort(centerAdjNodes));
                                 adjLeafNodes = curLeafNodes(adjacentNodes);
                                 validAdjacentNodes = cellfun(@(x) nnz(ismembc(x, centerAdjNodes)), adjLeafNodes) > 0;
                                 adjacentNodeNeigh(validAdjacentNodes) = min(adjacentNodeNeigh(validAdjacentNodes), neighItr);
                            end
                            
                            % If we've collected enough nodes, exit.
                            if nnz(validAdjacentNodes) == numel(validAdjacentNodes) || nnz(adjacentNodeNeigh == 255) > maxNodeDegree
  %                          if nnz(validAdjacentNodes) == numel(validAdjacentNodes)
                               break; 
                            end
                        end
                        validAdjacentNodes = adjacentNodeNeigh <= 3;
                        adjacentNodes = adjacentNodes(validAdjacentNodes);
                        adjacentNodeNeigh = adjacentNodeNeigh(validAdjacentNodes);
                    else
                        % For layer 1, we check for contour connectivity when
                        % forming an edge between two nodes.
                        centerRF = validMatrices{curNodeIds(nodeItr)};
                        adjacentNodeCoords = (curNodePreciseCoords(adjacentNodes,:) - centerPreciseArr(adjacentNodes,:)) + validMatrixHalfSize;
                        centerCoords = validMatrixHalfSize - (adjacentNodeCoords - validMatrixHalfSize); 
                        validAdjacentNodes = ones(numel(adjacentNodes),1) > 0;
                        shareRFCounts = zeros(numel(adjacentNodes),1);
                        for adjItr = 1:numel(adjacentNodes)
                            tempRF = centerRF(max(1, adjacentNodeCoords(adjItr,1) - incHalfMatrixSize):...
                                min(validMatrixSize, adjacentNodeCoords(adjItr,1) + incHalfMatrixSize), ...
                                max(1, adjacentNodeCoords(adjItr,2) - incHalfMatrixSize): ...
                                min(validMatrixSize, adjacentNodeCoords(adjItr,2) + incHalfMatrixSize));

                            % Now, we reverse the situation and get the part of
                            % the RF that is adjacent-node based.
                            adjRF = secValidMatrices{curNodeIds(adjacentNodes(adjItr))};
                            tempRF = tempRF & adjRF(max(1, centerCoords(adjItr,1) - incHalfMatrixSize):...
                                min(validMatrixSize, centerCoords(adjItr,1) + incHalfMatrixSize), ...
                                max(1, centerCoords(adjItr,2) - incHalfMatrixSize): ...
                                min(validMatrixSize, centerCoords(adjItr,2) + incHalfMatrixSize));

                            % If both rfs do not match, we move on.
                            if nnz(tempRF) == 0
                                validAdjacentNodes(adjItr) = 0;
                            else
                                shareRFCounts(adjItr) = nnz(tempRF);
                            end
                        end

                        % If this node has too many adjacent node pairs, we
                        % simply don't consider unlikely ones.
                        if nnz(shareRFCounts) > maxNodeDegree
                            [~, idx] = sort(shareRFCounts, 'descend');
                            adjacentNodes = adjacentNodes(sort(idx(1:maxNodeDegree)));
                        else
                            adjacentNodes = adjacentNodes(validAdjacentNodes);
                        end
                    end
                end
%                 
                if numel(adjacentNodes) == 0
                    imageEdgeDiagnostics(nodeItr) = 3;
                    continue;
                end
                
                %% Check for edge novelty.
                if currentLevelId > 1
                    % If edge novelty is enforced, act accordingly.
                    if edgeShareabilityThr < 1
                        commonLeafCounts = cellfun(@(x) sum(ismembc(x, curLeafNodes{nodeItr})), curLeafNodes(adjacentNodes));
                        novelNodes = (commonLeafCounts <= allowedSharedLeafNodes(adjacentNodes))' & ...
                            (commonLeafCounts <= allowedSharedLeafNodes(nodeItr))';
 %                       if nnz(novelNodes) > 1
                            adjacentNodes = adjacentNodes(novelNodes);
 %                           adjacentNodeNeigh = adjacentNodeNeigh(novelNodes);
 %                       else
 %                           novelNodes = (commonLeafCounts <= maxSharedLeafNodes(adjacentNodes))';
 %                           adjacentNodes = adjacentNodes(novelNodes);
 %                           adjacentNodeNeigh = adjacentNodeNeigh(novelNodes);
 %                       end
                    end
                end
                
                % Save info, if edge novelty threshold eliminated all
                % edges.
                if numel(adjacentNodes) == 0
                    imageEdgeDiagnostics(nodeItr) = 4;
                    continue;
                end
% 
%                 %% Now, we'll find a canonical edge representation for this node. 
%                 % We'll pick adjacent nodes one by one, based on their
%                 % contribution (addition of novel nodes), and their
%                 % distance. Ideally, we want to pick closer nodes with high
%                 % contributions.
%                 if currentLevelId > 1
%                      curNodeList = curLeafNodes{nodeItr};
%                      
%                      % Perform elimination of nodes that only has nodes
%                      % outside the RF.
% %                      if minimalEdgeCount
% %                          tempArr = curLeafNodes(adjacentNodes);
% %                          tempArr2 = cellfun(@(x) level1CoordsPooled(x,:), tempArr, 'UniformOutput', false);
% %                          tempArr = cellfun(@(x,y) x(y(:,1) > curNodeCoords(nodeItr,1) - rfRadius & ...
% %                               y(:,1) < curNodeCoords(nodeItr,1) + rfRadius & ...
% %                               y(:,2) > curNodeCoords(nodeItr,2) - rfRadius & ...
% %                               y(:,2) < curNodeCoords(nodeItr,2) + rfRadius), tempArr, tempArr2, 'UniformOutput', false);
% %                      else
%                           tempArr = curLeafNodes(adjacentNodes);
% %                     end
%                      
%                      if numel(adjacentNodes) > maxNodeDegree || minimalEdgeCount
% %                     if numel(adjacentNodes) > maxNodeDegree
%                          % Alternative 1.
%                          novelNodeCounts = cellfun(@(x) nnz(~ismembc(x, curNodeList)), tempArr)';
%                          if nnz(novelNodeCounts) > 0
%                              sortArr = [-novelNodeCounts, single(adjacentNodeNeigh), -distances(adjacentNodes),  single(curNodeIds(adjacentNodes))];
%                              [sortedArr, idx] = sortrows(sortArr);
%                              if nnz(~sortedArr(:,1)) == 0
%                                   limitVal = min([maxNodeDegree, size(sortedArr,1)]);
%                              else
%                                   limitVal = min([maxNodeDegree, size(sortedArr,1), find(~sortedArr(:,1), 1, 'first') - 1]);
%                              end
%                              adjacentNodes = adjacentNodes(sort(idx(1:limitVal)));
%                          else
%                              adjacentNodes = [];
%                          end
%                          
%                          if numel(adjacentNodes) == 0
%                               imageEdgeDiagnostics(nodeItr) = 5;
%                               continue;
%                          end
%                          
%                          % Alternative 2.
% %                            selectedNodeCount = 0;
% %                           pickedAdjacentNodes = [];
% %                           while selectedNodeCount < maxNodeDegree && ~isempty(adjacentNodes)
% %                               % Count the number of novel nodes for every adjacent
% %                               % node (which wasn't picked).
% %                               novelNodeCounts = cellfun(@(x) nnz(~ismembc(x, curNodeList)), curLeafNodes(adjacentNodes))';
% %                               if nnz(novelNodeCounts) == 0
% %                                    break;
% %                               end
% %                               
% %                               % Eliminate instances that have already been
% %                               % used for coverage.
% %                               validArr = novelNodeCounts > 0;
% % 
% %                               % Now, we pick next best adjacent node.
% %                       %        sortArr = [-novelNodeCounts, distances(adjacentNodes)];
% %                               sortArr = [-novelNodeCounts, curNodeIds(adjacentNodes)];
% %                               [~, idx] = sortrows(sortArr);
% %                               pickedAdjacentNodes = cat(1, pickedAdjacentNodes, adjacentNodes(idx(1)));
% %                               validArr(idx(1)) = 0;
% %                               curNodeList = fastsortedunique(sort(cat(2, curNodeList, curLeafNodes{adjacentNodes(idx(1))})));
% %                               selectedNodeCount = selectedNodeCount + 1;
% %                               
% %                               % Update lists.
% %                               adjacentNodes = adjacentNodes(validArr);
% %                           end
% %                           adjacentNodes = pickedAdjacentNodes;
%                      end
%                 else
%                    %% Eliminate far away nodes if this one has too many neighbors
%                    % Calculate scores (distances).
%                    scores = distances(adjacentNodes);
% 
%                    % Eliminate nodes having lower scores.
%                   if numel(adjacentNodes)>maxNodeDegree
%                         [idx] = getLargestNElements(scores, maxNodeDegree);
%                         adjacentNodes = adjacentNodes(idx);
%                   end
%                 end
                
                %% Assign final adjacent nodes.
                numberOfAdjacentNodes = numel(adjacentNodes);
                numberAdjArr(nodeItr) = numberOfAdjacentNodes;
                if numberOfAdjacentNodes > 0
                    curAdjacentNodes(nodeItr) = {adjacentNodes}; 
                else
                    curAdjacentNodes(nodeItr) = {double(emptyArr)};
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
                edgeDiagnostics(imageItr) = {imageEdgeDiagnostics};
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
            edgeDiagnostics(imageItr) = {imageEdgeDiagnostics};
            removedEdgeCounts = removedEdgeCounts + nnz(edgeIds == 0);
         end
         setGraphNodeSets{setItr} = imageGraphNodeSets;
         setEdgeDiagnostics{setItr} = edgeDiagnostics;
    end
    currentLevel = cat(1, setGraphNodeSets{:});
    currentLevel = [currentLevel{:}];
    
    % Obtain all diagnostics.
    edgeDiagnostics = cat(1, setEdgeDiagnostics{:});
    edgeDiagnostics = cat(1, edgeDiagnostics{:});
    
    %% Collect statistics about number of edges per node and print the stuff.
    edgeCounts = {currentLevel.adjInfo};
    edgeCounts = cellfun(@(x) size(x,1), edgeCounts);
    if usejava('jvm')
        figure('Visible', 'off'), hist(edgeCounts, 0:max(edgeCounts));
        saveas(gcf, [options.debugFolder '/level' num2str(currentLevelId) 'EdgHist.png']); 
    end
    display([num2str(nnz(edgeCounts == 0)) ' (%' num2str(100*nnz(edgeDiagnostics>0)/numel(edgeDiagnostics)) ') nodes have no edges.']);
    display([num2str(nnz(edgeDiagnostics == 1)) ' (%' num2str(100 * nnz(edgeDiagnostics == 1)/numel(edgeDiagnostics)) ') nodes are the only representatives in their images.']);
    display([num2str(nnz(edgeDiagnostics == 2)) ' (%' num2str(100 * nnz(edgeDiagnostics == 2)/numel(edgeDiagnostics)) ') nodes do not have edges cause there are no nodes in their RF.']);
    display([num2str(nnz(edgeDiagnostics == 3)) ' (%' num2str(100 * nnz(edgeDiagnostics == 3)/numel(edgeDiagnostics)) ') nodes do not have edges because of connectivity constraints.']);
    display([num2str(nnz(edgeDiagnostics == 4)) ' (%' num2str(100 * nnz(edgeDiagnostics == 4)/numel(edgeDiagnostics)) ') nodes do not have edges because of edge novelty constraint.']);
    display([num2str(nnz(edgeDiagnostics == 5)) ' (%' num2str(100 * nnz(edgeDiagnostics == 5)/numel(edgeDiagnostics)) ') nodes do not have edges because of enforcing RF limits in connections.']);
end
