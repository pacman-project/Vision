%> Name: learnModes
%>
%> Description: Given the node list of all images, this function learns the
%> modes by clustering pairwise relative positions in 2D space. 
%>
%> @param mainGraph The object graphs' data structure.
%> @param options Program options.
%> @param currentLevel The current scene graph level.
%> @param datasetName Name of the dataset.
%> 
%> @retval modes The mode list representing edge categories.
%>               modes are of the form: [ node11, node12, coord11, coord12;
%>                                        node21, node22, coord21, coord22;
%>                                      ...]
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 21.01.2014
function [modes] = learnModes(mainGraph, options, currentLevelId)
    display(['Learning modes for level ' num2str(currentLevelId) '...']);
    maxSamplesPerMode = options.mode.maxSamplesPerMode;
    minSamplesPerMode = options.mode.minSamplesPerMode;
    maximumModes = options.maximumModes;
    %% Step 0: Create initial data structures and initialize them.
    
    % Calculate edge radius.
    scale = (1/options.scaling)^(currentLevelId-1);
    neighborhood = fix(options.edgeRadius * scale);    
    
    % Set initial data structures for processing 
    currentLevel = mainGraph{currentLevelId};
    nodeIds = [currentLevel.labelId]';
    nodeCoords = cat(1, currentLevel.position);
    imageIds = [currentLevel.imageId]';
    
    %% Step 1: For each composition pair, get the 2-D distribution of samples.
    % Right now, we only work with 2-dimensional spatial relations.
    %% Put each image's node set into a different bin.
    numberOfImages = max(imageIds);
    imageGraphNodeSets = cell(numberOfImages, 1);
    imageNodeIdArr = cell(numberOfImages,1);
    imageNodeCoordArr = cell(numberOfImages,1);
    for imageItr = 1:max(imageIds)
       imageNodeIdx = imageIds == imageItr;
       imageGraphNodeSets(imageItr) = {currentLevel(imageIds == imageItr)};
       imageNodeIdArr(imageItr) = {nodeIds(imageNodeIdx)};
       imageNodeCoordArr(imageItr) = {nodeCoords(imageNodeIdx,:)};
    end
    
    %% Process each image separately (and in parallel)
    allSamples = cell(numberOfImages,1);
    parfor imageItr = 1:numberOfImages
        imageNodeIdx = imageIds == imageItr;

        % If there are no nodes in this image, move on.
        if nnz(imageNodeIdx) == 0
           continue; 
        end
        
        % Get data structures containing information about the nodes in this image.
        curNodeIds = imageNodeIdArr{imageItr};
        curNodeCoords = imageNodeCoordArr{imageItr};
        imageNodeIdx = find(imageNodeIdx)';
        numberOfNodes = numel(imageNodeIdx);
        
        %% Find all edges within this image.
        D = squareform(pdist(double(curNodeCoords)));
        D(logical(eye(size(D)))) = inf;
        D = D <= neighborhood;
        nodeStartIdx = 1;
        numberOfAllEdges = nnz(D);
        
        % If no edges exist, finish processing for this image.
        if numberOfAllEdges == 0
           continue;
        end
        
        % Assign the edges with the adjacency info.
        allEdges = zeros(nnz(D),2);
        for nodeItr = 1:numberOfNodes
           adjacentNodes = find(D(:,nodeItr));
           numberOfAdjacentNodes = numel(adjacentNodes);
           if numberOfAdjacentNodes==0
              continue; 
           end
           allEdges(nodeStartIdx:(nodeStartIdx+(numberOfAdjacentNodes-1)),:) = [repmat(nodeItr, numberOfAdjacentNodes,1), adjacentNodes];
           nodeStartIdx = nodeStartIdx + numberOfAdjacentNodes;
        end
        
        %% Fill in edge info With the labels of endpoints and their relative positions.
        node1Labels = curNodeIds(allEdges(:,1));
        node2Labels = curNodeIds(allEdges(:,2));
        node1Coords = curNodeCoords(allEdges(:,1),:);
        node2Coords = curNodeCoords(allEdges(:,2),:);
        edgeCoords = node1Coords - node2Coords;
        allSamples(imageItr) = {[node1Labels, node2Labels, edgeCoords]};
    end
    
    %% We have all possible edges extracted from all images. 
    % Moving on to mode calculation.
    allEdges = cat(1, allSamples{:});
    
    % If no edges exist, return.
    if isempty(allEdges)
       modes = [];
       return; 
    end

    % Eliminate edges for which second node's labels are smaller than first.
    allEdges = allEdges(allEdges(:,1) <= allEdges(:,2),:);
    
    %% Learn unique edge types and put them in cells for fast processing using parfor.
    [uniqueEdgeTypes, ~, IA] = unique(allEdges(:,1:2), 'rows');
    numberOfUniqueEdges = size(uniqueEdgeTypes,1);
    modes = cell(numberOfUniqueEdges,1);
    uniqueEdgeSamples = cell(numberOfUniqueEdges,1);
    parfor uniqueEdgeItr = 1:numberOfUniqueEdges
        samplesForEdge = allEdges(IA==uniqueEdgeItr,3:4); %#ok<PFBNS>
        %% If there are too many samples, get random samples.
        if size(samplesForEdge,1)>maxSamplesPerMode
            samplesForEdge = datasample(samplesForEdge, maxSamplesPerMode, 'Replace', false);
        end
        uniqueEdgeSamples(uniqueEdgeItr) = {samplesForEdge};
    end
    
    %% For each unique edge type (node1-node2 pair), estimate modes and save them in modes array.
    parfor uniqueEdgeItr = 1:numberOfUniqueEdges
  %      display(num2str(uniqueEdgeItr));
        w = warning('off', 'all');
        samples = double(uniqueEdgeSamples{uniqueEdgeItr});
        edgeType = uniqueEdgeTypes(uniqueEdgeItr,:);
        
        %% Assign a label to each sample.
        if maximumModes == 1
            classes = ones(size(samples,1),1);
        else
            classes = assignModes(samples, minSamplesPerMode, maximumModes);
        end

        %% Estimate cluster centers.
        numberOfClusters = max(classes);
        centers = zeros(numberOfClusters,4);
        centers(:,1:2) = repmat(edgeType, numberOfClusters, 1);
        
        for centerItr = 1:numberOfClusters
          clusterSamples = samples(classes==centerItr,:);
          centers(centerItr,3:4) = mean(clusterSamples,1);
        end
        
        modes(uniqueEdgeItr) = {centers};
           
        warning(w);
    end
    modes = round(cat(1, modes{:}));
        
    %% Add reverse modes to the modes array.
    reversedModes = modes(modes(:,1) ~= modes(:,2),:);
    tempArr = reversedModes(:,1);
    reversedModes(:,1) = reversedModes(:,2);
    reversedModes(:,2) = tempArr;
    reversedModes(:,3:4) = reversedModes(:,3:4) * -1;
    modes = [modes; reversedModes];

    % Sort array.
    modes = sortrows(modes);
    
end
