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
function [mainGraph] = createEdgesWithLabels(mainGraph, options, currentLevelId)
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
    matrixSize = [options.receptiveFieldSize, options.receptiveFieldSize];
    edgeType = options.edgeType;
    imagesPerSet = 100;
    
    %% Program options into variables.
    edgeNoveltyThr = 1-options.edgeNoveltyThr;
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
    
    %% Process each set separately (and in parallel)
    for setItr = 1:numberOfSets
         imageIdx = sets{setItr};
         imageNodeIdxSets = setNodeIdxSets{setItr};
         imageGraphNodeSets = setGraphNodeSets{setItr};
         imageNodeIdArr = setNodeIdArr{setItr} ;
         imageNodeCoordArr = setNodeCoordArr{setItr};
         disp(['Processing set ' num2str(setItr) ', which has ' num2str(numel(imageIdx)) ' images.']);
         
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

             %% Find all edges within this image.
             for nodeItr = 1:numberOfNodes
                centerArr = repmat(curNodeCoords(nodeItr,:), numberOfNodes,1);
                distances = zeros(numberOfNodes,1);
                adjacentNodes = curNodeCoords(:,1) > (centerArr(:,1) - halfMatrixSize) & ...
                     curNodeCoords(:,1) < (centerArr(:,1) + halfMatrixSize) & ...
                     curNodeCoords(:,2) > (centerArr(:,2) - halfMatrixSize) & ...
                     curNodeCoords(:,2) < (centerArr(:,2) + halfMatrixSize);
                distances(adjacentNodes) = sum((curNodeCoords(adjacentNodes,:) - centerArr(adjacentNodes,:)).^2, 2);

                %% Check for edge novelty.
                adjacentNodeIdx = find(adjacentNodes);
                centerLeafNodes = [];
                if currentLevelId > 1
                    centerLeafNodes = curLeafNodes{nodeItr};
                    commonLeafCounts = cellfun(@(x) sum(ismembc(x, centerLeafNodes)), curLeafNodes(adjacentNodes));
                    novelNodes = (commonLeafCounts <= maxSharedLeafNodes(adjacentNodes))';
                    adjacentNodes = adjacentNodeIdx(novelNodes);
                else
                    adjacentNodes = adjacentNodeIdx;
                end
                adjacentNodes = adjacentNodes(adjacentNodes~=nodeItr);

                %% A further check is done to ensure boundary/surface continuity.
                if strcmp(edgeType, 'continuity') && currentLevelId > 1
                    centerLeafEdges = cat(1, firstLevelAdjInfo{centerLeafNodes}); %#ok<PFBNS>
                    if isempty(centerLeafEdges)
                         adjacentNodes = [];
                    else
                         centerAdjNodes = sort(centerLeafEdges(:,2));
                         adjLeafNodes = curLeafNodes(adjacentNodes);
                         validAdjacentNodes = cellfun(@(x) nnz(ismembc(x, centerAdjNodes)), adjLeafNodes) > 0;
                         adjacentNodes = adjacentNodes(validAdjacentNodes);
                    end
                end           

                %% Now, we'll find a canonical edge representation for this node. 
                % We'll pick adjacent nodes one by one, based on their
                % contribution (addition of novel nodes), and their
                % distance. Ideally, we want to pick closer nodes with high
                % contributions.
                if currentLevelId > 1
                     curNodeList = centerLeafNodes;
                     unpickedAdjacentNodes = ones(numel(adjacentNodes),1) > 0;
                     selectedNodeCount = 0;
                     
                     while selectedNodeCount < maxNodeDegree
                         % If all nodes are picked, break.
                         if nnz(unpickedAdjacentNodes) == 0
                              break;
                         end

                         % Get the adjacent nodes to work with.
                         tempAdjacentNodes = adjacentNodes(unpickedAdjacentNodes);

                         % Count the number of novel nodes for every adjacent
                         % node (which wasn't picked).
                         novelNodeCounts = cellfun(@(x) nnz(~ismembc(x, curNodeList)), curLeafNodes(tempAdjacentNodes))';
                         if nnz(novelNodeCounts) == 0
                              break;
                         end

                         % Now, we pick next best adjacent node. 
                         sortArr = [-novelNodeCounts, distances(tempAdjacentNodes)];
                         [~, idx] = sortrows(sortArr);
                         unpickedAdjacentNodes(idx(1)) = 0;
                         curNodeList = fastsortedunique(sort(cat(2, curNodeList, curLeafNodes{tempAdjacentNodes(idx(1))})));
                         selectedNodeCount = selectedNodeCount + 1;
                     end
                     adjacentNodes = adjacentNodes(~unpickedAdjacentNodes);
                else
                      %% Eliminate far away nodes if this one has too many neighbors
                      % Calculate scores (distances).
                      scores = -distances(adjacentNodes);

                      % Eliminate nodes having lower scores.
                     if numel(adjacentNodes)>maxNodeDegree
                           [idx] = getLargestNElements(scores, maxNodeDegree);
                           adjacentNodes = adjacentNodes(idx);
                    end
                end

                %% Assign final adjacent nodes.
                curAdjacentNodes(nodeItr) = {[repmat(int32(nodeItr), numel(adjacentNodes),1), adjacentNodes]}; 
             end

             % Get rid of empty entries in curAdjacentNodes.
             nonemptyCurAdjacentNodeIdx = cellfun(@(x) ~isempty(x), curAdjacentNodes);
             curAdjacentNodes = curAdjacentNodes(nonemptyCurAdjacentNodeIdx);

             % Obtain edges and count them.
             allEdges = cat(1, curAdjacentNodes{:});
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
                for nodeItr = 1:numberOfNodes
                    edgeIdx = edges(:,1) == (nodeItr + imageNodeOffset);
                    if nnz(edgeIdx) > 0
                        curGraphNodes(nodeItr).adjInfo = edges(edgeIdx,:);
                    end
                end
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
