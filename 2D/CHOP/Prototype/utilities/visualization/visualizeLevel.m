%> Name: visualizeLevel
%>
%> Description: Given the current and previous vocabulary level, this
%> function visualizes the current vocabulary level as a separate patche for
%> each word in the vocabulary level.
%>
%> @param currentLevel Current vocabulary level.
%> @param graphLevel Current graph level.
%> @param levelId Identifier of the current level.
%> @param numberOfPrevNodes Number of nodes in previous vocabulary level.
%> @param options Program options.
%> @param isRedundant If currentLevel consists of redundant compositions,
%> set 1. Otherwise set 0.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 10.02.2014
%> Redundant vocabulary output option added. 10.05.2014
function [allNodeInstances, representativeNodes] = visualizeLevel( currentLevel, vocabulary, graphLevel, firstActivations, leafNodes, leafNodeCoords, levelId, numberOfFirstLevelNodes, ~, options)
    % Read options to use in this file.
    currentFolder = options.currentFolder;
    datasetName = options.datasetName;
    filtBandCount = size(options.filters{1},3);
    numberOfThreads = options.numberOfThreads;
    childrenPerNode = options.vis.nodeReconstructionChildren;
    instancePerNode = options.vis.instancePerNode;
    visualizedNodes = options.vis.visualizedNodes;
    vocabLevelLabels = [currentLevel.label];
    
    % Process filters for visualization, remove zero values.
    filters = options.filters;
    filters = cellfun(@(x) (x - min(min(x))) / (max(max(x)) - min(min(x))), filters, 'UniformOutput', false);
    for filterItr = 1:numel(filters)
         curFilter = filters{filterItr};
         curFilter(curFilter<0.001) = 0.001;
         filters{filterItr} = curFilter;
    end
    
    vocabLevelIdx = zeros(max(vocabLevelLabels),1);
    for itr = 1:max(vocabLevelLabels)
         vocabLevelIdx(itr) = find(vocabLevelLabels == itr, 1, 'first');
    end
    allNodeInstances = [];
    multiNodeInstances = {currentLevel.adjInfo};
    multiNodeInstances = cellfun(@(x) ~isempty(x), multiNodeInstances);
    
    % For the first level, we only print a single instance.
    if levelId == 1
        instanceImgDim = 1; 
    else
        instanceImgDim = round(sqrt(instancePerNode));
    end
    
    % We decrease this number by 1, since the best match is printed twice.
    instancePerNode = instancePerNode-1;
    filterType = options.filterType;
    if strcmp(filterType, 'auto') && levelId == 1
        isAutoFilter = true;
        deadFeatures = options.auto.deadFeatures;
    else
        isAutoFilter = false;
    end
    
    %% Learn positions of leaf nodes.
    if levelId > 1
        leafNodeSets = {graphLevel.leafNodes}';
        centerPos = int32(round(cat(1, graphLevel.precisePosition)));
        nodeLabelIds = double(cat(1, graphLevel.realLabelId));
        orNodeLabelIds = double(cat(1, graphLevel.labelId));
        nodeImageIds = cat(1, graphLevel.imageId);
        nodeActivations = cat(1, graphLevel.activation);
        leafNodeLabelIds = leafNodes(:, 1);
        leafNodePos = leafNodeCoords;
    end
    
    % Find representative parts for every or node (the one with most
    % instances).
    if levelId > 1
         representativeNodes = zeros(max(vocabLevelLabels),1);
         orNodeLabelIds = double(cat(1, graphLevel.labelId));
         for nodeItr = 1:numel(representativeNodes)
              representativeNodes(nodeItr) = mode(nodeLabelIds(orNodeLabelIds == nodeItr));
              if ~multiNodeInstances(representativeNodes(nodeItr))
                   instances = nodeLabelIds(orNodeLabelIds == nodeItr);
                   instances = instances(multiNodeInstances(instances));
                   if ~isempty(instances)
                        representativeNodes(nodeItr) = mode(instances);
                   end
              end
         end
    else
         representativeNodes = (1:numel(currentLevel))';
    end
    
    %% Create output folder structure.
    reconstructionDir = [currentFolder '/debug/' datasetName '/level' num2str(levelId) '/reconstruction/'];
    if ~exist(reconstructionDir, 'dir')
       mkdir(reconstructionDir);
    else
       rmdir(reconstructionDir, 's');
       mkdir(reconstructionDir);
    end
    
    %% In level 1, only print low level filters as first n nodes.
    if levelId == 1
        filterDir = [currentFolder '/filters/' options.filterType '/'];
        numberOfNodes = double(max(vocabLevelLabels));
        nodeImgs = cell(numberOfNodes,1);
        for nodeItr = 1:numberOfNodes
            mask = imread([filterDir 'filt' num2str(nodeItr) '.png']);
            nodeImgs{nodeItr} = {mask};
            imwrite(mask, [reconstructionDir num2str(nodeItr) '.png']);
        end
    else
        %% In other levels, combine the nodes of previous levels depending on mode info and visualize current level.
        % Read previous layer's masks.
        firstLevelDir = [currentFolder '/debug/' datasetName '/level' num2str(1) '/reconstruction/'];
%        numberOfNodes = double(max(vocabLevelLabels));
        numberOfNodes = numel(vocabLevelLabels);
        firstNodeMasks = cell(numberOfFirstLevelNodes,1);
        avgFirstNodeMasks = cell(numberOfFirstLevelNodes,1);
        firstLevelPatchLowDims = zeros(numberOfFirstLevelNodes,2);
        firstLevelPatchHighDims = zeros(numberOfFirstLevelNodes,2);
        
        %% Collect patches for the first level (leaf nodes).
        for nodeItr = 1:numberOfFirstLevelNodes
            tempImg = double(imread([firstLevelDir num2str(nodeItr) '.png']));
            tempImg = tempImg / 255;
            firstNodeMasks(nodeItr) = {tempImg};
            avgFirstNodeMasks(nodeItr) = {mean(tempImg,3)};
            firstLevelPatchHighDims(nodeItr,:) = floor([size(firstNodeMasks{nodeItr},1), size(firstNodeMasks{nodeItr},2)]/2);
            firstLevelPatchLowDims(nodeItr,:) = ([size(firstNodeMasks{nodeItr},1), size(firstNodeMasks{nodeItr},2)] - firstLevelPatchHighDims(nodeItr,:)) - 1;
        end
        
        %% To parallelize things, we put vocabulary nodes in different sets, and give each to a thread.
        nodeSet = (1:numberOfNodes);
        numberOfThreadsUsed = 1;
        if options.parallelProcessing && numberOfThreads > 1 && numberOfNodes > numberOfThreads
            parallelNodeSetIdx = round(1:((numberOfNodes-1)/numberOfThreads):numberOfNodes);
            parallelNodeSetIdx = parallelNodeSetIdx(2:end)-parallelNodeSetIdx(1:(end-1));
            parallelNodeSetIdx(end) = parallelNodeSetIdx(end) + 1;
            parallelNodeSets = mat2cell(nodeSet, 1, parallelNodeSetIdx);
            parallelVocabNodeSets = parallelNodeSets;
            numberOfThreadsUsed = numberOfThreads;
        else
            parallelNodeSets = {nodeSet};
            parallelVocabNodeSets = {1:numberOfNodes};
        end
        
        %% Go through each node in the current layer, and reconstuct it to
        % get its mask in the end. Each node is reconstructed using the
        % nodes in the previous layer which contribute to its definition. 
        setImgs = cell(numberOfThreadsUsed,1);
        allNodeInstances = cell(numberOfNodes,1);
        for setItr = 1:numberOfThreadsUsed
            w = warning('off', 'all');
            nodeSet = parallelNodeSets{setItr};
            vocabNodeSet = parallelVocabNodeSets{setItr};
            nodeImgs = cell(numel(nodeSet),1);
            
            % Go through each composition in current node set.
            for nodeItr = 1:numel(nodeSet)
                %% Get the children (leaf nodes) from all possible instance in the dataset. Keep the info.
   %             labelId = vocabLevelIdx(vocabNodeSet(nodeItr));
                labelId = vocabNodeSet(nodeItr);
          
                orNodeId = find(representativeNodes == labelId);
                if ismember(labelId, representativeNodes)
                     nodeInstances = find(orNodeLabelIds==orNodeId);
                     validIdx = orNodeLabelIds==orNodeId;
                else
                     nodeInstances = find(nodeLabelIds == labelId);
                     validIdx = nodeLabelIds == labelId;
                end
                instanceActivations = nodeActivations(validIdx);
                instanceImageIds = nodeImageIds(validIdx);
                bestNodeInstance = -1;
                if ~isempty(instanceActivations)
                     [~, sortIdx] = sort(instanceActivations, 'descend');
                     if numel(nodeInstances)>instancePerNode
                          % Get best instances to print.
                          [~, validSortIdx] = unique(instanceImageIds(sortIdx), 'stable');
                          sortIdx = sortIdx(validSortIdx);
                          if numel(sortIdx) > instancePerNode
                              sortIdx = sortIdx(1:instancePerNode);
                          end
                          nodeInstances = sort(nodeInstances(sortIdx));
                     end
                end
                nodeInstances = [bestNodeInstance; nodeInstances]; %#ok<AGROW>
                instanceImgs = cell(numel(nodeInstances),1);
                nodeInstances = [nodeInstances, zeros(numel(nodeInstances),4)]; %#ok<AGROW>
                
                %% Now, we print the rest of the instances.
                for nodeInstanceItr = 1:size(nodeInstances,1)
                     if nodeInstanceItr == 1
                         % It is supposed to provide an approximate view that the algorithm has learned.
                         currentMask = imread([options.debugFolder '/level' num2str(levelId) '/modalProjection/' num2str(labelId) '.png']);
                         imgSize = [size(currentMask,1), size(currentMask,2)];
                     else
                         nodeInstance = nodeInstances(nodeInstanceItr,1);
                         instanceLeafNodes = leafNodeSets{nodeInstance};
                         projectedNodes = zeros(numel(instanceLeafNodes), 4);
                         projectedNodes(:,1) = double(leafNodeLabelIds(instanceLeafNodes));
                         projectedNodes(:,2:3) = double(leafNodeCoords(instanceLeafNodes,:));
                         projectedNodes(:,4) = firstActivations(instanceLeafNodes,:);
                         
                         % Update coordinates.
                         mins = min(projectedNodes(:,2:3), [], 1);
                         maxs = max(projectedNodes(:,2:3), [], 1);
                         meanCoordOffset = (mins + maxs) / 2;
                         coordOffset = repmat(imgSize/2  - meanCoordOffset, size(projectedNodes,1), 1);
                         projectedNodes(:,2:3) = round(projectedNodes(:,2:3) + coordOffset);
                         bounds = [1,1,imgSize] - [round(coordOffset(1,1:2)), round(coordOffset(1,1:2))];
                         nodeInstances(nodeInstanceItr,2:end) = bounds;
                         
                         % Obtain a PoE visualization
                         [currentMask, ~, ~, ~] = obtainPoE(projectedNodes, [], [], [], imgSize, filters);
                         currentMask = uint8(currentMask*255);
                     end
                     
                    %% Print the files to output folders.
                    if nodeInstanceItr == 1
                        % Create a false color image.
                        imwrite(currentMask, [reconstructionDir num2str(nodeSet(nodeItr)) '.png']);
                    else
                        imwrite(currentMask, [reconstructionDir num2str(nodeSet(nodeItr)) '_var_' num2str(nodeInstanceItr) '.png']);
                    end
                    
                    % Save all instances. We'll print them to another
                    % image.
                    instanceImgs(nodeInstanceItr) = {currentMask};
                end
                
                % Combine instance images together, write them to a larger
                % image and save it.
                allNodeInstances{labelId} = nodeInstances;
                nodeImgs(nodeItr) = {instanceImgs};
            end
            warning(w);
            setImgs(setItr) = {nodeImgs};
        end
        nodeImgs = cat(1, setImgs{:});
        clear leafNodeSets centerPos nodeLabelIds leafNodeLabelIds leafNodePos prevNodeMasks avgPrevNodeMasks parallelNodeSets;
    end
    
    %% Combine all compositions and show them within a single image.
    if levelId > 1
         visRepresentativeNodes = representativeNodes(multiNodeInstances(representativeNodes));
    else
         visRepresentativeNodes = representativeNodes;
    end
    % Learn number of rows/columns, use only representative nodes.
    numberOfNodes = numel(visRepresentativeNodes);
    colImgCount = ceil(sqrt(numberOfNodes));
    rowImgCount = ceil(numberOfNodes/colImgCount);
    if numberOfNodes < visualizedNodes
       visualizedNodes = numberOfNodes; 
    end
    smallColImgCount = ceil(sqrt(visualizedNodes));
    smallRowImgCount = ceil(visualizedNodes/smallColImgCount);

    %% Create large images that have multiple components.
    dim3 = size(nodeImgs{1}{1},3);
    imgSize = [size(nodeImgs{1}{1},1), size(nodeImgs{1}{1},2)];
    overallImage = NaN((rowImgCount)*(imgSize(1)+1)+1, colImgCount * (imgSize(2)+1)+1, dim3);
    overallInstanceImage = NaN((smallRowImgCount * instanceImgDim)*(imgSize(1)+1)+1, smallColImgCount * instanceImgDim * (imgSize(2)+1)+1, dim3);
 %   overallInstanceRealImage = NaN((rowImgCount * instanceImgDim)*(compMaskSize(1)+1)+1, colImgCount * instanceImgDim * (compMaskSize(2)+1)+1, dim3);
    for nodeItr = 1:numberOfNodes
        instanceImgs = nodeImgs{visRepresentativeNodes(nodeItr)};
        for instItr = 1:numel(instanceImgs)
            compFinalMask = instanceImgs{instItr};

            % If this feature is a dead one, reduce illumination by three.
            if isAutoFilter && ismember(visRepresentativeNodes(nodeItr), deadFeatures)
                compFinalMask = uint8(round(compFinalMask * 0.33));
            end

            % Add the composition's mask to the overall mask image.
            rowStart2= 2 + floor((nodeItr-1)/smallColImgCount)*(imgSize(1)+1) *instanceImgDim ;
            colStart2 = 2 + rem(nodeItr-1, smallColImgCount) * (imgSize(2)+1) * instanceImgDim;
            if instItr == 1
                rowStart = 2 + floor((nodeItr-1)/colImgCount)*(imgSize(1)+1);
                colStart = 2 + rem(nodeItr-1, colImgCount) * (imgSize(2)+1);

                overallImage(rowStart:(rowStart+imgSize(1)-1), ...
                    colStart:(colStart+imgSize(2)-1), :) = instanceImgs{1};
            end
            %We're writing sample instances to other images. Find where to
            %write them and put them in their location.
            rowInstStart = floor((instItr - 1)/instanceImgDim)*(imgSize(1)+1);
            colInstStart = rem(instItr-1, instanceImgDim) * (imgSize(2)+1);
            if nodeItr<=visualizedNodes
                overallInstanceImage((rowStart2 + rowInstStart):((rowStart2 + rowInstStart)+imgSize(1)-1), ...
                    (colStart2+colInstStart):((colStart2+colInstStart)+imgSize(2)-1), :) = compFinalMask;
            end
        end
    end

    clear instanceImgs nodeImgs setImgs;

    % A final make up in order to separate masks from each other by 1s.
    whiteRowIdx = 1:(imgSize(1)+1) *instanceImgDim:size(overallInstanceImage,1);
    whiteColIdx = 1:(imgSize(2)+1) *instanceImgDim:size(overallInstanceImage,2);
    overallImage(isnan(overallImage)) = 255;
    overallImage = uint8(overallImage);
    overallInstanceImage(isnan(overallInstanceImage)) = 0;
    overallInstanceImage(:, whiteColIdx, :) = 255;
    overallInstanceImage(whiteRowIdx, :, :) = 255;
    overallInstanceImage = uint8(overallInstanceImage);

    % Then, write the compositions the final image.
    imwrite(overallImage, [currentFolder '/debug/' datasetName '/level' num2str(levelId) '_vb.png']);
    imwrite(overallInstanceImage, [currentFolder '/debug/' datasetName '/level' num2str(levelId) '_vb_variations.png']);
end



