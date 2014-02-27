%> Name: applyLocalInhibition
%>
%> Description: Applies local inhibition of nodes to eliminate those who do
%> not introduce novelty, into the graph. In this context, novelty is
%> measured by the ratio of leaf nodes not existing in the main graph over
%> all leaf nodes of the tested node. In case of an overlap, the node whose
%> label id is better wins. graphLevel is ASSUMED to be ordered by
%> labelIds.
%>
%> @param vocabLevel The current vocabulary level, ASSUMED ordered by
%> mdlValues.
%> @param graphLevel The current graph level, ASSUMED ordered by labelIds.
%> @param currentModes Current modes representing geometric relationships.
%> @param options Program options.
%> @param levelItr Current level number.
%>
%> @retval graphLevel Modified graph level.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 11.12.2013
function [graphLevel] = applyLocalInhibition(vocabLevel, graphLevel, currentModes, options, levelItr)
    % Calculate edge radius.
    scale = (1/options.scaling)^(levelItr-1);
    neighborhood = fix(options.edgeRadius * scale);
    noveltyThr = options.noveltyThr;

    %% First inhibition is done via structural evaluation. 
    % If the structures of two different subs match, it is very likely that
    % they have same node sets. We must eliminate such cases.
    % We only need mode positions.
    if options.useReceptiveField
        currentModes = currentModes(:,3:4);
        currentModes((size(currentModes,1)+1),:) = [0, 0];

        % Get unique sub descriptions independent of the center node. Each
        % description is a string that encodes coarse structure of the sub.
        vocabNodeLabels = {vocabLevel.children};
        vocabEdges = {vocabLevel.adjInfo};
        newMode = size(currentModes,1);
        vocabNeighborModes = cellfun(@(x) [newMode; x(:,3)], vocabEdges, 'UniformOutput', false);
        vocabNodePositions = cellfun(@(x) currentModes(x,:) - repmat(min(currentModes(x,:)), numel(x), 1), vocabNeighborModes, 'UniformOutput', false);
        
        % Get rid of noise in node positions to have a coarse geometric
        % representation.
        normalizeConstant = 4/max(max(currentModes));
        vocabNodePositions = cellfun(@(x) round(x*normalizeConstant), vocabNodePositions, 'UniformOutput', false);
        
        vocabSortOrder = cell(size(vocabLevel,1),1);
        for vocabNodeItr = 1:numel(vocabLevel)
            [~, vocabSortOrder{vocabNodeItr}] = sortrows(vocabNodePositions{vocabNodeItr});
        end
        vocabDescriptions = cellfun(@(x,y,z) [x(z)', y(z,:)], vocabNodeLabels, vocabNodePositions, vocabSortOrder, 'UniformOutput', false);
        vocabDescriptions = cellfun(@(x) mat2str(x(:)'), vocabDescriptions, 'UniformOutput', false);
    
        % Now, we have unique vocabulary descriptions, independent of center
        % nodes. Eliminate those match, by keeping only one (the one with
        % highest mdlScore) of each set.
 %       vocabDescriptions = vocabDescriptions(numel(vocabLevel):-1:1);
        [~, IA, ~] = unique(vocabDescriptions, 'stable');
        labelIds = cat(1, graphLevel.labelId);
        graphLevel = graphLevel(ismember(labelIds, IA));
    end

    % Fill in necessary internal structures.
    imageIds = [graphLevel.imageId];
    numberOfImages = max(imageIds);
    
    % Get node coordinates.
    nodeCoords = cat(1, graphLevel.position);
    
    %% Prepare data structures for parallel processing.
    imageGraphLevels = cell(numberOfImages,1);
    imageAllNodeCoords = cell(numberOfImages,1);
    for imageId = 1:numberOfImages
        imageNodeIdx = imageIds == imageId;
        imageGraphLevels(imageId) = {graphLevel(1,imageNodeIdx)};
        imageAllNodeCoords(imageId) = {nodeCoords(imageNodeIdx,:)}; 
    end
    
    %% Go over each node and check neighboring nodes for novelty introduced. Eliminate weak ones.
    preservedNodes = cell(numberOfImages,1);
    parfor imageId = 1:numberOfImages
        imageGraphLevel = imageGraphLevels{imageId};
        imageNodeCoords = imageAllNodeCoords{imageId};
        numberOfNodesInImage = numel(imageGraphLevel);
        imagePreservedNodes = ones(numberOfNodesInImage,1);
        imageLeafNodes = {imageGraphLevel.leafNodes}';
        maxSharedLeafNodes = cellfun(@(x) numel(x) * noveltyThr , imageLeafNodes, 'UniformOutput', false);
        
        for nodeItr = 1:(numberOfNodesInImage-1)
          %% If nobody has erased this node before, it has a right to be in the final graph.
          if imagePreservedNodes(nodeItr) == 0
              continue;
          end
          
          %% Get each neighboring node.
          thisNodeCoords = imageNodeCoords(nodeItr,:);
          centerArr = repmat(thisNodeCoords, numberOfNodesInImage, 1);
          distances = sqrt(sum((centerArr - imageNodeCoords).^2, 2));
          adjacentNodes = distances <= neighborhood; 
          adjacentNodes(1:nodeItr) = 0;
          selfLeafNodes = imageLeafNodes{nodeItr};
          
          %% Go over each adjacent node, and apply inhibition if their leaf nodes are too common under current novelty threshold.
          imagePreservedNodes(adjacentNodes) = cellfun(@(x,y) sum(ismembc(x, selfLeafNodes)) <= y, ...
              imageLeafNodes(adjacentNodes), maxSharedLeafNodes(adjacentNodes));
        end
        preservedNodes(imageId) = {imagePreservedNodes'};
    end
    preservedNodes = [preservedNodes{:}]>0;
    graphLevel = graphLevel(:,preservedNodes);
end
