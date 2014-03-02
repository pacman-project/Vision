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
    useReceptiveField = options.useReceptiveField;
    maxSamplesPerMode = options.mode.maxSamplesPerMode;
    %% Step 0: Create initial data structures and initialize them.
    
    % Calculate edge radius.
    scale = (1/options.scaling)^(currentLevelId-1);
    neighborhood = fix(options.edgeRadius * scale);
    
    % Set initial data structures for processing 
    currentLevel = mainGraph{currentLevelId};
    nodeIds = [currentLevel.labelId]';
    nodeCoords = cat(1, currentLevel.position);
    imageIds = [currentLevel.imageId]';
    numberOfNodes = max(nodeIds);
    
    %% Step 1: For each composition pair, get the 2-D distribution of samples.
    % Right now, we only work with 2-dimensional spatial relations.
    
    % Get node image ids and coords into separate arrays to reduce overhead
    % in parfor.
    nodeImageIdArr = cell(numberOfNodes,1);
    nodeCoordArr = cell(numberOfNodes,1);
    for nodeItr = 1:numberOfNodes
       nodeIdx = nodeIds==nodeItr;
       nodeImageIdArr(nodeItr) = {imageIds(nodeIdx, :)};
       nodeCoordArr(nodeItr) = {nodeCoords(nodeIdx, :)};
    end
    
    modes = cell(numberOfNodes,1);
    parfor node1 = 1:numberOfNodes
       w = warning('off', 'all');
       firstNodeImageIds = nodeImageIdArr{node1};
       firstNodeCoords = nodeCoordArr{node1};
       curModes = cell(numberOfNodes,1);
       
       %    Prevent empty cluster warnings in kmeans.
       for node2 = node1:numberOfNodes
           % Get first type of nodes.
           
           % Get second type of node ids and positions.
           secondNodeImageIds = nodeImageIdArr{node2};
           secondNodeCoords = nodeCoordArr{node2};
           
           % Get common image ids.
           commonImageIds = fastintersect(fastsortedunique(firstNodeImageIds), ...
               fastsortedunique(secondNodeImageIds));
           
           % If no common imageids exist, nothing to do here.
           if isempty(commonImageIds)
              continue; 
           end
           
           % Find all edges in between.
           samples = cell(numel(commonImageIds),1);
           for imageIdItr = 1:numel(commonImageIds)
               % Get center coords belonging to first node type in image.
               imageId = commonImageIds(imageIdItr);
               firstNodeCenterCoords = firstNodeCoords(firstNodeImageIds == imageId,:);
               adjacentSecondNodeCoords = secondNodeCoords(secondNodeImageIds == imageId,:);
               numberOfCenters = size(firstNodeCenterCoords,1);
               numberOfAdjacentNodes = size(adjacentSecondNodeCoords,1);
               
               % Replicate second node coords as many times as number of
               % centers.
               firstNodeCenterCoords = kron(firstNodeCenterCoords,ones(numberOfAdjacentNodes,1));
               adjacentSecondNodeCoords = repmat(adjacentSecondNodeCoords, numberOfCenters,1);
               
               % Get relative coordinates of both sets of nodes.
               samples(imageIdItr) = {adjacentSecondNodeCoords - firstNodeCenterCoords};
           end
           
           % Combine all samples, and get rid of invalid ones.
           samples = cat(1,samples{:});
           distances = sqrt(sum(samples.^2,2));
           if node1==node2
              samples = samples(distances > 0 & distances <= neighborhood,:);
              
              % If receptive fields are not used, we only use ~half of the
              % samples on the right side of a separating line, if node1 == node2.
              sampleSums = sum(samples,2);
              if ~useReceptiveField
                  samples = samples((samples(:,1) >= 0 & sampleSums >= 0) | ...
                      (samples(:,1) < 0 & sampleSums > 0), :);
              end
           else
              samples = samples(distances <= neighborhood,:);
           end
           
           % If no samples remain, continue.
           if isempty(samples)
              continue; 
           end
           
           %% If there are too many samples, get random samples.
           if size(samples,1)>maxSamplesPerMode
                samples = datasample(samples, maxSamplesPerMode, 'Replace', false);
           end
           
           %% Assign a label to each sample.
           classes = assignModes(samples, options);
           
           %% Estimate cluster centers.
           numberOfClusters = max(classes);
           centers = zeros(numberOfClusters,4);
           centers(:,1:2) = repmat([node1, node2], numberOfClusters, 1);
           for centerItr = 1:numberOfClusters
              centers(centerItr,3:4) = mean(samples(classes==centerItr,:),1);
           end
           
           % Assign centers to modes.
           curModes(node2) = {centers};

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
                   imwrite(label2rgb(distributionImg, 'jet', 'k', 'shuffle'), ...
                       [options.currentFolder '/debug/' datasetName '/level' num2str(currentLevelId) '/pairwise/' num2str(node1) '_' num2str(node2) '.png']);
               end
           end
       end
       warning(w);
       modes(node1) = {curModes};
    end
    %% If modes have been found, concatenate and return them.
    modes = cat(1, modes{:});
    modes = fix(cat(1, modes{:}));
    
    %% If receptive fields are used, add reverse modes to the modes array.
    if options.useReceptiveField
        reversedModes = modes(modes(:,1) ~= modes(:,2),:);
        tempArr = reversedModes(:,1);
        reversedModes(:,1) = reversedModes(:,2);
        reversedModes(:,2) = tempArr;
        reversedModes(:,3:4) = reversedModes(:,3:4) * -1;
        modes = [modes; reversedModes];
        
        % Sort array.
        modes = sortrows(modes);
    end
end