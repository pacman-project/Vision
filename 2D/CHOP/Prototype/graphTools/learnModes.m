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
function [modes] = learnModes(mainGraph, options, currentLevelId, datasetName)
    display(['Learning modes for level ' num2str(currentLevelId) '...']);
    %% Step 0: Create initial data structures and initialize them.
    % Prevent empty cluster warnings in kmeans.
    w = warning('off', 'all');
    
    % Calculate edge radius.
    scale = (1/options.scaling)^(currentLevelId-1);
    neighborhood = fix(options.edgeRadius * scale);
    
    % Set initial data structures for processing 
    currentLevel = mainGraph{currentLevelId};
    nodeIds = [currentLevel.labelId]';
    nodeCoords = cat(1, currentLevel.position);
    imageIds = [currentLevel.imageId]';
    rfIds = [currentLevel.rfId]';
    isCenter = [currentLevel.isCenter]';
    numberOfNodes = max(nodeIds);
    imageRfPairs = [imageIds, rfIds];
    
    %% Step 1: For each composition pair, get the 2-D distribution of samples.
    % Right now, we only work with 2-dimensional spatial relations.
    modeOffset = 0;
    modeIdxOffset = 1;
    modes = cell(numberOfNodes * numberOfNodes,1);
    
    for node1 = 1:numberOfNodes
       for node2 = node1:numberOfNodes
           % Get first type of nodes.
           firstNodeIdx = find(nodeIds==node1);
           firstNodeImageRfPairs = imageRfPairs(firstNodeIdx, :);
           
           % Get second type of nodes.
           secondNodeIdx = find(nodeIds==node2);
           secondNodeImageRfPairs = imageRfPairs(secondNodeIdx, :);
           
           % Get common image - rf pairs.
           commonImageRfPairs = intersect(firstNodeImageRfPairs, ...
               secondNodeImageRfPairs, 'rows');
           
           % If no common image-rf pairs exist, nothing to do here.
           if isempty(commonImageRfPairs)
              continue; 
           end
           
           firstNodeCoords = nodeCoords(firstNodeIdx, :);
           firstNodeIsCenter = isCenter(firstNodeIdx,:);
           secondNodeCoords = nodeCoords(secondNodeIdx, :);
           
           % Find all edges in between.
           samples = cell(size(commonImageRfPairs,1),1);
           for imageRfPairId = 1:size(commonImageRfPairs,1)
               % Get center coords belonging to first node type in image.
               imageRfPair = commonImageRfPairs(imageRfPairId,:);
               firstNodeCenterCoords = firstNodeCoords(firstNodeImageRfPairs(:,1) == imageRfPair(1) & firstNodeImageRfPairs(:,2) == imageRfPair(2) & firstNodeIsCenter,:);
               adjacentSecondNodeCoords = secondNodeCoords(secondNodeImageRfPairs(:,1) == imageRfPair(1) & secondNodeImageRfPairs(:,2) == imageRfPair(2),:);
               numberOfCenters = size(firstNodeCenterCoords,1);
               numberOfAdjacentNodes = size(adjacentSecondNodeCoords,1);
               
               % Replicate second node coords as many times as number of
               % centers.
               firstNodeCenterCoords = kron(firstNodeCenterCoords,ones(numberOfAdjacentNodes,1));
               adjacentSecondNodeCoords = repmat(adjacentSecondNodeCoords, numberOfCenters,1);
               
               % Get relative coordinates of both sets of nodes.
               samples(imageRfPairId) = {adjacentSecondNodeCoords - firstNodeCenterCoords};
           end
           
           % Combine all samples, and get rid of invalid ones.
           samples = cat(1,samples{:});
           distances = sqrt(sum(samples.^2,2));
           if node1==node2
              samples = samples(distances > 0 & distances <= neighborhood,:);
              
              % If receptive fields are not used, we only use ~half of the
              % samples on the right side of a separating line, if node1 == node2.
              sumSamples = sum(samples,2);
              if ~options.useReceptiveField
                  samples = samples((samples(:,1) >= 0 & sumSamples >= 0) | ...
                      (samples(:,1) < 0 & sumSamples > 0), :);
              end
           else
              samples = samples(distances <= neighborhood,:);
           end
           
           %% Assign a label to each sample.
           classes = assignModes(samples, options);
           
           %% Estimate cluster centers.
           numberOfClusters = max(classes);
           centers = zeros(numberOfClusters,4);
           centers(:,1) = ones(numberOfClusters,1) * node1;
           centers(:,2) = ones(numberOfClusters,1) * node2;
           for centerItr = 1:numberOfClusters
              centers(centerItr,3:4) = mean(samples(classes==centerItr,:),1);
           end
           
           % Assign centers to modes.
           modes(modeIdxOffset) = {centers};
           modeIdxOffset = modeIdxOffset + 1;

           % Change mode offset and move on.
           modeOffset = modeOffset + max(classes);
           
           %% In debug mode, write classes to the output as images.
           if options.debug
               distributionImg = zeros(options.maxImageDim*2+1);
               samplesToWrite = floor(samples + options.maxImageDim + 1);
               
               % If no samples are to be written, move on.
               if numel(samplesToWrite) < 1
                   continue;
               end
                   
               samplesInd = sub2ind(size(distributionImg), samplesToWrite(:,1), samplesToWrite(:,2));
               distributionImg(samplesInd) = classes;

               % Resize the distribution image so it is of the smallest
               % possible size.
               bound = max(max(abs(samples)));
               midPoint = options.maxImageDim + 1;
               distributionImg = distributionImg((midPoint-bound):(midPoint+bound), (midPoint-bound):(midPoint+bound));
               if ~exist([options.currentFolder '/debug/' datasetName '/level' num2str(currentLevelId) '/pairwise/'], 'dir')
                   mkdir([options.currentFolder '/debug/' datasetName '/level' num2str(currentLevelId) '/pairwise/']);
               end
               if ~isempty(distributionImg)
                   imwrite(label2rgb(distributionImg, 'jet', 'k', 'shuffle'), [options.currentFolder '/debug/' datasetName '/level' num2str(currentLevelId) '/pairwise/' num2str(node1) '_' num2str(node2) '.png']);
               end
           end
       end
    end
    %% If modes have been found, concatenate and return them.
    if modeIdxOffset > 1
        modes = cat(1, modes{1:(modeIdxOffset-1),:});
    else
        modes = 0;
    end
    warning(w);
end