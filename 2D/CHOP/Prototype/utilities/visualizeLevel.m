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
function [] = visualizeLevel( currentLevel, vocabulary, graphLevel, firstActivations, leafNodes, levelId, numberOfFirstLevelNodes, numberOfPrevNodes, options)
    % Read options to use in this file.
    currentFolder = options.currentFolder;
    datasetName = options.datasetName;
    filtBandCount = size(options.filters{1},3);
    numberOfThreads = options.numberOfThreads;
    childrenPerNode = options.vis.nodeReconstructionChildren;
    instancePerNode = options.vis.instancePerNode;
    visualizedNodes = options.vis.visualizedNodes;
    
    % These parameter relate to drawing the approximate model of each node
    % in the vocabulary.
    edgeIdMatrix = options.edgeIdMatrix;
    midPointEdgeIdMatrix = int32(round(size(edgeIdMatrix) / 2));
    edgeIdMatrixSize = size(edgeIdMatrix);
    scale = (1/options.scaling)^(levelId-2);
    neighborhood = floor(options.edgeRadius * scale);
    upsampleRatio = neighborhood / floor((options.edgeQuantize-1)/2);
    if strcmp(options.filterType, 'gabor')
       inhibitionRadius = options.gabor.inhibitionRadius; 
    else
       inhibitionRadius = options.auto.inhibitionRadius;
    end
    inhibitionRadiusSq = inhibitionRadius^2;
    
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
        centerPos = cat(1, graphLevel.position);
        nodeLabelIds = cat(1, graphLevel.labelId);
        nodeImageIds = cat(1, graphLevel.imageId);
        nodeActivations = cat(1, graphLevel.activation);
        leafNodeLabelIds = cat(1, leafNodes(:, 1));
        leafNodePos = cat(1,leafNodes(:,2:3));
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
        numberOfNodes = numel(currentLevel);
        nodeImgs = cell(numberOfNodes,1);
        for nodeItr = 1:numberOfNodes
            mask = double(imread([filterDir 'filt' num2str(nodeItr) '.png']));
            mask = (mask - min(min(min(mask)))) / (max(max(max(mask))) - min(min(min(mask))));
            nodeImgs{nodeItr} = {uint8(round(mask * 255))};
            imwrite(mask, [reconstructionDir num2str(nodeItr) '.png']);
        end
    else
        %% In other levels, combine the nodes of previous levels depending on mode info and visualize current level.
        % Read previous layer's masks.
        firstLevelDir = [currentFolder '/debug/' datasetName '/level' num2str(1) '/reconstruction/'];
        prevLevelDir = [currentFolder '/debug/' datasetName '/level' num2str(levelId-1) '/reconstruction/'];
        numberOfNodes = numel(currentLevel);
        firstNodeMasks = cell(numberOfFirstLevelNodes,1);
        prevNodeMasks = cell(numberOfPrevNodes,1);
        avgFirstNodeMasks = cell(numberOfFirstLevelNodes,1);
        avgPrevNodeMasks = cell(numberOfPrevNodes,1);
        firstLevelPatchLowDims = zeros(numberOfFirstLevelNodes,2);
        firstLevelPatchHighDims = zeros(numberOfFirstLevelNodes,2);
        prevLevelPatchLowDims = zeros(numberOfPrevNodes,2);
        prevLevelPatchHighDims = zeros(numberOfPrevNodes,2);
        
        %% Collect patches for the first level (leaf nodes).
        for nodeItr = 1:numberOfFirstLevelNodes
            tempImg = double(imread([firstLevelDir num2str(nodeItr) '.png']));
            tempImg = (tempImg - min(min(min(tempImg)))) / (max(max(max(tempImg))) - min(min(min(tempImg))));
            firstNodeMasks(nodeItr) = {tempImg};
            avgFirstNodeMasks(nodeItr) = {mean(tempImg,3)};
            firstLevelPatchHighDims(nodeItr,:) = floor([size(firstNodeMasks{nodeItr},1), size(firstNodeMasks{nodeItr},2)]/2);
            firstLevelPatchLowDims(nodeItr,:) = ([size(firstNodeMasks{nodeItr},1), size(firstNodeMasks{nodeItr},2)] - firstLevelPatchHighDims(nodeItr,:)) - 1;
%            lowResponseThrs(nodeItr) = max(max(max(tempImg)))/10;
        end
        
        %% Collect patches for the previous level (for approximate model printing).
        if levelId > 2
             for nodeItr = 1:numberOfPrevNodes
                   tempImg = double(imread([prevLevelDir num2str(nodeItr) '.png']));
                   tempImg = (tempImg - min(min(min(tempImg)))) / (max(max(max(tempImg))) - min(min(min(tempImg))));
                   prevNodeMasks(nodeItr) = {tempImg};
                   avgPrevNodeMasks(nodeItr) = {mean(tempImg,3)};
                   prevLevelPatchHighDims(nodeItr,:) = floor([size(prevNodeMasks{nodeItr},1), size(prevNodeMasks{nodeItr},2)]/2);
                   prevLevelPatchLowDims(nodeItr,:) = ([size(prevNodeMasks{nodeItr},1), size(prevNodeMasks{nodeItr},2)] - prevLevelPatchHighDims(nodeItr,:)) - 1;
             end
        else
             prevNodeMasks = firstNodeMasks;
             avgPrevNodeMasks = avgFirstNodeMasks;
             prevLevelPatchHighDims = firstLevelPatchHighDims;
             prevLevelPatchLowDims = firstLevelPatchLowDims;
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
        for setItr = 1:numberOfThreadsUsed
            w = warning('off', 'all');
            nodeSet = parallelNodeSets{setItr};
            vocabNodeSet = parallelVocabNodeSets{setItr};
            nodeImgs = cell(numel(nodeSet),1);
            
            % Go through each composition in current node set.
            for nodeItr = 1:numel(nodeSet)
                %% Get the children (leaf nodes) from all possible instance in the dataset. Keep the info.
                labelId = vocabNodeSet(nodeItr);
                nodeInstances = find(nodeLabelIds==labelId);
                validIdx = nodeLabelIds==labelId;
                instanceActivations = nodeActivations(validIdx);
                instanceImageIds = nodeImageIds(validIdx);
                [~, sortIdx] = sort(instanceActivations, 'descend');
                bestNodeInstance = nodeInstances(sortIdx(1));
                if numel(nodeInstances)>instancePerNode
                     % Get best instances to print.
                     [~, validSortIdx] = unique(instanceImageIds(sortIdx), 'stable');
                     sortIdx = sortIdx(validSortIdx);
                     if numel(sortIdx) > instancePerNode
                         sortIdx = sortIdx(1:instancePerNode);
                     end
                     nodeInstances = sort(nodeInstances(sortIdx));
                end
                nodeInstances = [bestNodeInstance; nodeInstances]; %#ok<AGROW>
                
                instanceImgs = cell(numel(nodeInstances),1);
                
                %% Now, we print the rest of the instances.
                for nodeInstanceItr = 1:numel(nodeInstances)
                     if nodeInstanceItr == 1
                         %% Here, we get the description of the node, and print that.
                                                  %% Here, we get the description of the node, and print that.
%                              % It is supposed to provide an approximate view that the algorithm has learned.
%                             currentLevel(labelId).children;
%                             children = (currentLevel(labelId).children)';
%                             childrenCoords = zeros(numel(children), 2, 'int32');
%                             % Find the location of every child based on the
%                             % first one.
%                             nodeAdjInfo = currentLevel(labelId).adjInfo;
%                             for childItr = 2:numel(children)
%                                  [childrenCoords(childItr, 1), childrenCoords(childItr,2)] = ind2sub(edgeIdMatrixSize, find(edgeIdMatrix == nodeAdjInfo(childItr-1, 3)));
%                                  childrenCoords(childItr,:) = childrenCoords(childItr,:) - midPointEdgeIdMatrix;
%                             end
%                             % We scale the coordinates by a factor based on the
%                             % current level id. Then, we round them to get 
%                             % actual image coordinates (relative still).
%                             childrenCoords = int32(round(double(childrenCoords) * upsampleRatio));
% 
%                             %% Then, we perform a simple inhibition process to remove overlapping level 1 instances.
%                             validNodes = ones(size(children,1),1) > 0;
%                             numberOfAllNodes = size(children,1);
%                             for nodeItr2 = 1:(numberOfAllNodes-1)
%                                 if ~validNodes(nodeItr2) 
%                                    continue; 
%                                 end
%                                 remainingCoords = childrenCoords((nodeItr2+1):end, :);
%                                 validRemainingNodes = sum((remainingCoords - repmat(childrenCoords(nodeItr2,:), ...
%                                     numberOfAllNodes - nodeItr2, 1)).^2, 2) >= inhibitionRadiusSq;
%                                 validNodes((nodeItr2+1):end) = validNodes((nodeItr2+1):end) & validRemainingNodes;
%                             end
%                             children = children(validNodes);
%                             childrenCoords = childrenCoords(validNodes,:);
%     
%                             %% Finally, we use previous level's data structures for printing.
%                             patchLowDims = prevLevelPatchLowDims;
%                             patchHighDims = prevLevelPatchHighDims;
%                             avgNodeMasks = avgPrevNodeMasks;
%                             nodeMasks = prevNodeMasks;
                         
                         % It is supposed to provide an approximate view that the algorithm has learned.
                         projectedNodes = projectNode([labelId, 0, 0, levelId], vocabulary, inhibitionRadius);
                         children = projectedNodes(:,1);
                         childrenCoords = projectedNodes(:,2:3);
                         patchLowDims = firstLevelPatchLowDims;
                         patchHighDims = firstLevelPatchHighDims;
                         avgNodeMasks = avgFirstNodeMasks;
                         nodeMasks = firstNodeMasks;
                     else
                         nodeInstance = nodeInstances(nodeInstanceItr);
                         instanceLeafNodes = leafNodeSets{nodeInstance};
                         instancePos = mat2cell(centerPos(nodeInstance,:), ones(1, numel(nodeInstance)), 2);
                         instanceLeafNodeSets = leafNodeSets(nodeInstance,:);
                         instanceLeafNodePos = cellfun(@(x, y) leafNodePos(x,:) - ...
                             repmat(y, numel(x), 1), instanceLeafNodeSets, instancePos, ...
                             'UniformOutput', false);
                         children = leafNodeLabelIds(cat(2, instanceLeafNodeSets{:}));
                         childrenCoords= cat(1, instanceLeafNodePos{:});

                         % Trim children if total number is more than a threshold.
                         if numel(children) > childrenPerNode
                            children = children(1:childrenPerNode,:);
                            childrenCoords = childrenCoords(1:childrenPerNode,:);
                         end
                     end
                     
                     %% We use first level's data structures for printing, when printing the instances.
                     if nodeInstanceItr == 2
                        patchLowDims = firstLevelPatchLowDims;
                        patchHighDims = firstLevelPatchHighDims;
                        avgNodeMasks = avgFirstNodeMasks;
                        nodeMasks = firstNodeMasks;
                     end
                     
                    %% At this point, we have the relative coordinates of all children. 
                    % All we will do is to create an empty mask large enough, and
                    % write the children's masks over it.
                    % Determine the upper, lower, left and right bounds (extremes)
                    % of the mask required to hold this composition's mask.
                    maskMinX = 0;
                    maskMinY = 0;
                    maskMaxX = 0;
                    maskMaxY = 0;
                    for childItr = 1:numel(children)
                        minX = childrenCoords(childItr,1) - patchLowDims(children(childItr),1);
                        maxX = childrenCoords(childItr,1) + patchHighDims(children(childItr),1);
                        minY = childrenCoords(childItr,2) - patchLowDims(children(childItr),2);
                        maxY = childrenCoords(childItr,2) + patchHighDims(children(childItr),2);
                        if maskMinX > minX
                            maskMinX = minX;
                        end
                        if maskMaxX < maxX
                            maskMaxX = maxX;
                        end
                        if maskMinY > minY
                            maskMinY = minY;
                        end
                        if maskMaxY < maxY
                            maskMaxY = maxY;
                        end
                    end
                    maskMinX = double(floor(maskMinX)-1);
                    maskMinY = double(floor(maskMinY)-1);
                    maskMaxX = double(ceil(maskMaxX)+1);
                    maskMaxY = double(ceil(maskMaxY)+1);
                    childrenCoords = int32(round(double(childrenCoords) - [ones(numel(children),1) * maskMinX, ones(numel(children),1) * maskMinY]));

                    %% Write the children's masks to the current mask.
                    currentMask = zeros((maskMaxX - maskMinX)+1, (maskMaxY - maskMinY)+1, filtBandCount);
                    currentFilledMask = zeros((maskMaxX - maskMinX)+1, (maskMaxY - maskMinY)+1, filtBandCount);
                    currentLabelMask = zeros((maskMaxX - maskMinX)+1, (maskMaxY - maskMinY)+1);
                    for childItr = 1:numel(children)
                        % Write the child's mask to the output.
                        if nodeInstanceItr == 1
                            assignedWeight = 1;
                        else
                            assignedWeight = firstActivations(instanceLeafNodes(childItr));
                        end
                        
                        currentMask((childrenCoords(childItr,1)-patchLowDims(children(childItr),1)):(childrenCoords(childItr,1)+patchHighDims(children(childItr),1)), ...
                          (childrenCoords(childItr,2)-patchLowDims(children(childItr),2)):(childrenCoords(childItr,2)+patchHighDims(children(childItr),2)),:) = ...
                          max(assignedWeight * nodeMasks{children(childItr)}, ...
                              (currentMask((childrenCoords(childItr,1)-patchLowDims(children(childItr),1)):(childrenCoords(childItr,1)+patchHighDims(children(childItr),1)), ...
                                  (childrenCoords(childItr,2)-patchLowDims(children(childItr),2)):(childrenCoords(childItr,2)+patchHighDims(children(childItr),2)),:)));

                        currentFilledMask((childrenCoords(childItr,1)-patchLowDims(children(childItr),1)):(childrenCoords(childItr,1)+patchHighDims(children(childItr),1)), ...
                              (childrenCoords(childItr,2)-patchLowDims(children(childItr),2)):(childrenCoords(childItr,2)+patchHighDims(children(childItr),2))) = 1;

                        % Mark label image for the child.
                        currentLabelMask((childrenCoords(childItr,1)-patchLowDims(children(childItr),1)):(childrenCoords(childItr,1)+patchHighDims(children(childItr),1)), ...
                             (childrenCoords(childItr,2)-patchLowDims(children(childItr),2)):(childrenCoords(childItr,2)+patchHighDims(children(childItr),2))) = ...
                                     max(currentLabelMask((childrenCoords(childItr,1)-patchLowDims(children(childItr),1)):(childrenCoords(childItr,1)+patchHighDims(children(childItr),1)), ...
                             (childrenCoords(childItr,2)-patchLowDims(children(childItr),2)):(childrenCoords(childItr,2)+patchHighDims(children(childItr),2))), ...
                             double(avgNodeMasks{children(childItr)}>0.2) * childItr);
                    end

                    %% Add background to currentMask, and normalize it.
                    % Learn the median color to use as background.
                    validValues = currentMask(currentFilledMask>0);
                    filledValue = median(validValues);

                    % Assign filling value to each band.
                    for bandItr = 1:size(currentMask,3)
                        bandMask = currentMask(:,:,bandItr);
                        bandMask(~currentLabelMask) = filledValue;
                        currentMask(:,:,bandItr) = bandMask;
                    end

                    % Normalize current mask to [0,1].
                    currentMask = ( currentMask - min(validValues)) / (max(validValues) - min(validValues));

                    %% If the size of the current mask can be divided by two, pad sides to prevent this.
                    dimRems = rem(size(currentMask),2);
                    if dimRems(1) == 0
                        currentMask = [currentMask; ones(1, size(currentMask,2), size(currentMask,3)) * filledValue];
                        currentLabelMask = [currentLabelMask; ones(1, size(currentLabelMask,2), size(currentLabelMask,3)) * filledValue];
                    end
                    if dimRems(2) == 0
                        currentMask = [currentMask, ones(size(currentMask, 1), 1, size(currentMask,3)) * filledValue];
                        currentLabelMask = [currentLabelMask, ones(size(currentLabelMask, 1), 1, size(currentLabelMask,3)) * filledValue];
                    end

                    % Assign currentMask to the label image.
                    currentLabelImg = currentMask;
                    
                    %% Print the files to output folders.
                    if nodeInstanceItr == 1
                        % Create a false color image.
                        falseColorImg = double(label2rgb(int32(currentLabelMask), 'hot', 'k', 'shuffle'));
                        if isequal(size(falseColorImg), size(currentMask))
                             falseColorImg = uint8(round(falseColorImg .* currentMask));
                        else
                             for bandItr = 1:size(falseColorImg,3)
                                  falseColorImg(:,:,bandItr) = currentMask(:,:,1) .* falseColorImg(:,:,bandItr);
                             end
                             falseColorImg = uint8(round(falseColorImg));
                        end
                        
                        imwrite(currentMask, [reconstructionDir num2str(nodeSet(nodeItr)) '.png']);
                        imwrite(falseColorImg, [reconstructionDir num2str(nodeSet(nodeItr)) '_falseColor.png']);
                        imwrite(currentMask, [reconstructionDir num2str(nodeSet(nodeItr)) '_' num2str(currentLevel(nodeSet(nodeItr)).mdlScore) '.png']);
                        imwrite(currentLabelImg, [reconstructionDir num2str(nodeSet(nodeItr)) '_comp.png']);
                    else
                        imwrite(currentMask, [reconstructionDir num2str(nodeSet(nodeItr)) '_var_' num2str(nodeInstanceItr) '.png']);
                    end
                    
                    % Save all instances. We'll print them to another
                    % image.
                    currentLabelImg = uint8(round(currentLabelImg * 255));
                    instanceImgs(nodeInstanceItr) = {currentLabelImg};
                end
                % Combine instance images together, write them to a larger
                % image and save it.
                nodeImgs(nodeItr) = {instanceImgs};
            end
            warning(w);
            setImgs(setItr) = {nodeImgs};
        end
        nodeImgs = cat(1, setImgs{:});
        clear leafNodeSets centerPos nodeLabelIds leafNodeLabelIds leafNodePos prevNodeMasks avgPrevNodeMasks parallelNodeSets;
    end
    
    %% Combine all compositions and show them within a single image.
    % Learn number of rows/columns.
    colImgCount = ceil(sqrt(numberOfNodes));
    rowImgCount = ceil(numberOfNodes/colImgCount);
    if numberOfNodes < visualizedNodes
       visualizedNodes = numberOfNodes; 
    end
    smallColImgCount = ceil(sqrt(visualizedNodes));
    smallRowImgCount = ceil(visualizedNodes/smallColImgCount);

    %% Show the set of compositions in a single image and save it.
    % Read all masks into an array, and get the extreme dimensions.allCompMasks = cell(numberOfNodes,1);
    compMaskSize = [1, 1];
    for nodeItr = 1:numberOfNodes
        instanceImgs = nodeImgs{nodeItr};
        if ~isempty(instanceImgs)
            dim3 = size(instanceImgs{1},3);
        end
        instanceImgSizes = cellfun(@(x) [size(x,1), size(x,2)], instanceImgs, 'UniformOutput', false);
        if ~isempty(instanceImgSizes)
            compMaskSize = max(compMaskSize, max(cat(1, instanceImgSizes{:})));
        end
    end

    % Make mask sizes uniform and write them all back.
    if levelId>1
        for nodeItr = 1:numberOfNodes
            instanceImgs = nodeImgs{nodeItr};
            for instItr = 1:numel(instanceImgs)
                tempMask2 = instanceImgs{instItr};
                finalTempMask = zeros([compMaskSize, size(tempMask2,3)], 'uint8');
                margins = (compMaskSize - [size(tempMask2,1), size(tempMask2,2)])/2;
                finalTempMask((floor(margins(1))+1):(end-ceil(margins(1))), ...
                    (floor(margins(2))+1):(end-ceil(margins(2))), :) = tempMask2;
                % A make-up to fill in NaNs (empty points).
      %          fillInValue = 0;
                fillInValues = double(finalTempMask(finalTempMask>0 & finalTempMask<255));
                fillInValue = fillInValues(1);
           %     fillInValue = median(double(finalTempMask(finalTempMask>0 & finalTempMask<255)));
                finalTempMask(finalTempMask == 0) = fillInValue;
                instanceImgs{instItr} = finalTempMask;
            end
            nodeImgs{nodeItr} = instanceImgs;
        end
    end

    % Using the maximum dimensions, transform each composition image to the
    % same size. 
    overallImage = NaN((rowImgCount)*(compMaskSize(1)+1)+1, colImgCount * (compMaskSize(2)+1)+1, dim3);
    overallInstanceImage = NaN((smallRowImgCount * instanceImgDim)*(compMaskSize(1)+1)+1, smallColImgCount * instanceImgDim * (compMaskSize(2)+1)+1, dim3);
 %   overallInstanceRealImage = NaN((rowImgCount * instanceImgDim)*(compMaskSize(1)+1)+1, colImgCount * instanceImgDim * (compMaskSize(2)+1)+1, dim3);
    for nodeItr = 1:numberOfNodes
        instanceImgs = nodeImgs{nodeItr};
        for instItr = 1:numel(instanceImgs)
            compFinalMask = instanceImgs{instItr};

            % If this feature is a dead one, reduce illumination by three.
            if isAutoFilter && ismember(nodeItr, deadFeatures)
                compFinalMask = uint8(round(compFinalMask * 0.33));
            end

            % Add the composition's mask to the overall mask image.
            rowStart = 2 + floor((nodeItr-1)/colImgCount)*(compMaskSize(1)+1);
            colStart = 2 + rem(nodeItr-1, colImgCount) * (compMaskSize(2)+1);
            rowStart2= 2 + floor((nodeItr-1)/smallColImgCount)*(compMaskSize(1)+1) *instanceImgDim ;
            colStart2 = 2 + rem(nodeItr-1, smallColImgCount) * (compMaskSize(2)+1) * instanceImgDim;
            if instItr == 1
                overallImage(rowStart:(rowStart+compMaskSize(1)-1), ...
                    colStart:(colStart+compMaskSize(2)-1), :) = compFinalMask;
            end
            %We're writing sample instances to other images. Find where to
            %write them and put them in their location.
            rowInstStart = floor((instItr - 1)/instanceImgDim)*(compMaskSize(1)+1);
            colInstStart = rem(instItr-1, instanceImgDim) * (compMaskSize(2)+1);
            if nodeItr<=visualizedNodes
                overallInstanceImage((rowStart2 + rowInstStart):((rowStart2 + rowInstStart)+compMaskSize(1)-1), ...
                    (colStart2+colInstStart):((colStart2+colInstStart)+compMaskSize(2)-1), :) = compFinalMask;
            end
        end
        if ~isempty(instanceImgs)
            imwrite(instanceImgs{1}, [reconstructionDir num2str(nodeItr) '_uni.png']);
        end
    end

    clear instanceImgs nodeImgs setImgs;

    % A final make up in order to separate masks from each other by 1s.
    whiteRowIdx = 1:(compMaskSize(1)+1) *instanceImgDim:size(overallInstanceImage,1);
    whiteColIdx = 1:(compMaskSize(2)+1) *instanceImgDim:size(overallInstanceImage,2);
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



