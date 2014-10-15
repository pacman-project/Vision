%> Name: visualizeLevel
%>
%> Description: Given the current and previous vocabulary level, this
%> function visualizes the current vocabulary level as a separate patche for
%> each word in the vocabulary level.
%>
%> @param currentLevel Current vocabulary level.
%> @param graphLevel Current graph level.
%> @param levelId Identifier of the current level.
%> @param modes Modes of the previous level to reconstruct the features.
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
function [] = visualizeLevel( currentLevel, graphLevel, leafNodes, leafDistanceMatrix, levelId, ~, numberOfPrevNodes, options, isRedundant)
    % Read options to use in this file.
    currentFolder = options.currentFolder;
    datasetName = options.datasetName;
    filtBandCount = size(options.filters{1},3);
    numberOfThreads = options.numberOfThreads;
    childrenPerNode = options.vis.nodeReconstructionChildren;
    filterType = options.filterType;
    if strcmp(filterType, 'auto') && levelId == 1
        isAutoFilter = true;
        deadFeatures = options.auto.deadFeatures;
    else
        isAutoFilter = false;
    end
    
    % Turn leaf distance matrix into leaf similarity matrix.
    leafSimilarityMatrix = (max(max(leafDistanceMatrix)) - leafDistanceMatrix) / max(max(leafDistanceMatrix));
    
    %% Learn positions of leaf nodes.
    if levelId > 1
        leafNodeSets = {graphLevel.leafNodes}';
        centerPos = cat(1, graphLevel.position);
        nodeLabelIds = cat(1, graphLevel.labelId);
        leafNodeLabelIds = cat(1, leafNodes{:, 1});
        leafNodePos = cat(1,leafNodes{:,2});
    end
    
    %% Create output folder structure.
    reconstructionDir = [currentFolder '/debug/' datasetName '/level' num2str(levelId) '/reconstruction/'];
    if ~exist(reconstructionDir, 'dir')
       mkdir(reconstructionDir);
    else
       rmdir(reconstructionDir, 's');
       mkdir(reconstructionDir);
    end
    
    %% Read label ids if redundant level is processed.
    if isRedundant
        labelIds = [currentLevel.label];
    else
        labelIds = [];
    end
    
    %% In level 1, only print low level filters as first n nodes.
    if levelId == 1
        filterDir = [currentFolder '/filters/' options.filterType '/'];
        numberOfNodes = numel(currentLevel);
        for nodeItr = 1:numberOfNodes
            mask = double(imread([filterDir 'filt' num2str(nodeItr) '.png']));
            mask = (mask - min(min(min(mask)))) / (max(max(max(mask))) - min(min(min(mask))));
            imwrite(mask, [reconstructionDir num2str(nodeItr) '.png']);
        end
    else
        %% In other levels, combine the nodes of previous levels depending on mode info and visualize current level.
        % Read previous layer's masks.
        prevLevelDir = [currentFolder '/debug/' datasetName '/level' num2str(1) '/reconstruction/'];
        numberOfNodes = numel(currentLevel);
        prevNodeMasks = cell(numberOfPrevNodes,1);
        avgPrevNodeMasks = cell(numberOfPrevNodes,1);
        patchLowDims = zeros(numberOfPrevNodes,2);
        patchHighDims = zeros(numberOfPrevNodes,2);
        lowResponseThrs = zeros(numberOfPrevNodes,1);
        
        for nodeItr = 1:numberOfPrevNodes
            tempImg = double(imread([prevLevelDir num2str(nodeItr) '.png']));
            tempImg = (tempImg - min(min(min(tempImg)))) / (max(max(max(tempImg))) - min(min(min(tempImg))));
            prevNodeMasks(nodeItr) = {tempImg};
            avgPrevNodeMasks(nodeItr) = {mean(tempImg,3)};
            patchHighDims(nodeItr,:) = floor([size(prevNodeMasks{nodeItr},1), size(prevNodeMasks{nodeItr},2)]/2);
            patchLowDims(nodeItr,:) = ([size(prevNodeMasks{nodeItr},1), size(prevNodeMasks{nodeItr},2)] - patchHighDims(nodeItr,:)) - 1;
            lowResponseThrs(nodeItr) = max(max(max(tempImg)))/10;
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
%       for setItr = 1:numberOfThreadsUsed
        parfor setItr = 1:numberOfThreadsUsed
            w = warning('off', 'all');
            nodeSet = parallelNodeSets{setItr};
            vocabNodeSet = parallelVocabNodeSets{setItr};
            
            % Go through each composition in current node set.
            for nodeItr = 1:numel(nodeSet)
                %% Get the children (leaf nodes) from all possible instance in the dataset. Keep the info.
                labelId = vocabNodeSet(nodeItr);
                nodeInstances = find(nodeLabelIds==labelId);
                if numel(nodeInstances)>5
		     nodeInstances = nodeInstances(1:5);   % CHANGE: Print only the first realization.
		end
                instancePos = mat2cell(centerPos(nodeInstances,:), ones(1, numel(nodeInstances)), 2);
                instanceLeafNodeSets = leafNodeSets(nodeInstances,:);
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

                %% Give weight the children to get a cleaner representation of the average image for the node.
                % First, we learn base scores, which essentially describe
                % how close each child is to each other. A child does
                % contribute to itself.
                distMatrix = squareform(pdist(childrenCoords));
                maxDist = max(max(distMatrix));
                scoreMatrix = exp(maxDist - distMatrix);
                maxDist = max(max(scoreMatrix));
                scoreMatrix = scoreMatrix / maxDist;

                % Now, weight the scores by the similarity matrix entries.
                weightMatrix = leafSimilarityMatrix(children, children);

                % Multiply score matrix by weight matrix to get weighted
                % scores (element by element).
                scoreMatrix = scoreMatrix .* weightMatrix;
                scoreMatrix = sum(scoreMatrix,1);
                scoreMatrix = scoreMatrix / max(scoreMatrix);

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
                maskMinX = floor(maskMinX)-1;
                maskMinY = floor(maskMinY)-1;
                maskMaxX = ceil(maskMaxX)+1;
                maskMaxY = ceil(maskMaxY)+1;
                childrenCoords = round(childrenCoords - [ones(numel(children),1) * maskMinX, ones(numel(children),1) * maskMinY]);

                %% Write the children's masks to the current mask.
                currentMask = zeros((maskMaxX - maskMinX)+1, (maskMaxY - maskMinY)+1, filtBandCount);
                currentFilledMask = zeros((maskMaxX - maskMinX)+1, (maskMaxY - maskMinY)+1, filtBandCount);
                currentLabelMask = zeros((maskMaxX - maskMinX)+1, (maskMaxY - maskMinY)+1);
                for childItr = 1:numel(children)
                    % Write the child's mask to the output.
                    currentMask((childrenCoords(childItr,1)-patchLowDims(children(childItr),1)):(childrenCoords(childItr,1)+patchHighDims(children(childItr),1)), ...
                      (childrenCoords(childItr,2)-patchLowDims(children(childItr),2)):(childrenCoords(childItr,2)+patchHighDims(children(childItr),2)),:) = ...
                      max(prevNodeMasks{children(childItr)}.*scoreMatrix(childItr), ...
                          currentMask((childrenCoords(childItr,1)-patchLowDims(children(childItr),1)):(childrenCoords(childItr,1)+patchHighDims(children(childItr),1)), ...
                              (childrenCoords(childItr,2)-patchLowDims(children(childItr),2)):(childrenCoords(childItr,2)+patchHighDims(children(childItr),2)),:));

                    currentFilledMask((childrenCoords(childItr,1)-patchLowDims(children(childItr),1)):(childrenCoords(childItr,1)+patchHighDims(children(childItr),1)), ...
                          (childrenCoords(childItr,2)-patchLowDims(children(childItr),2)):(childrenCoords(childItr,2)+patchHighDims(children(childItr),2))) = 1;

                    % Mark label image for the child.
                    currentLabelMask((childrenCoords(childItr,1)-patchLowDims(children(childItr),1)):(childrenCoords(childItr,1)+patchHighDims(children(childItr),1)), ...
                         (childrenCoords(childItr,2)-patchLowDims(children(childItr),2)):(childrenCoords(childItr,2)+patchHighDims(children(childItr),2))) = ...
                                 max(currentLabelMask((childrenCoords(childItr,1)-patchLowDims(children(childItr),1)):(childrenCoords(childItr,1)+patchHighDims(children(childItr),1)), ...
                         (childrenCoords(childItr,2)-patchLowDims(children(childItr),2)):(childrenCoords(childItr,2)+patchHighDims(children(childItr),2))), ...
                         double(avgPrevNodeMasks{children(childItr)}>lowResponseThrs(children(childItr))) * childItr);
                end
                
                %% Add background to currentMask, and normalize it.
                % Learn the median color to use as background..
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
                end
                if dimRems(2) == 0
                    currentMask = [currentMask, ones(size(currentMask, 1), 1, size(currentMask,3)) * filledValue];
                end

                % Assign currentMask to the label image.
                currentLabelImg = currentMask;

                %% Print the files to output folders.
                if isRedundant
                    realLabel = labelIds(nodeSet(nodeItr));
                    optionOrder = nnz(labelIds(1:nodeSet(nodeItr))==realLabel);
                    imwrite(currentMask, [reconstructionDir num2str(realLabel) '_option' num2str(optionOrder) '.png']);
                else
                    imwrite(currentMask, [reconstructionDir num2str(nodeSet(nodeItr)) '.png']);
                    imwrite(currentLabelImg, [reconstructionDir num2str(nodeSet(nodeItr)) '_comp.png']);
                end
            end
            warning(w);
        end
    end
    
    %% Combine all compositions and show them wit`hin a single image.
    if ~isRedundant
        % Learn number of rows/columns.
        colImgCount = ceil(sqrt(numberOfNodes));
        rowImgCount = ceil(numberOfNodes/colImgCount);

        %% Show the set of compositions in a single image and save it.
        % Read all masks into an array, and get the extreme dimensions.
        allCompMasks = cell(numberOfNodes,1);
        compMaskSize = [1, 1];
        for nodeItr = 1:numberOfNodes
            % Read image and add it to the figure.
            if levelId == 1
                tempMask = imread([reconstructionDir num2str(nodeItr) '.png']);
            else
                tempMask = imread([reconstructionDir num2str(nodeItr) '_comp.png']);
            end
            compMaskSize = max(compMaskSize, [size(tempMask,1), size(tempMask,2)]);
            allCompMasks(nodeItr) = {tempMask};
        end
        
        % Make mask sizes uniform and write them all back.
        if levelId>1
            for nodeItr = 1:numberOfNodes
                tempMask2 = imread([reconstructionDir num2str(nodeItr) '.png']);
                finalTempMask = zeros([compMaskSize, size(tempMask2,3)], 'uint8');
                margins = (compMaskSize - [size(tempMask2,1), size(tempMask2,2)])/2;
                finalTempMask((floor(margins(1))+1):(end-ceil(margins(1))), ...
                    (floor(margins(2))+1):(end-ceil(margins(2))), :) = tempMask2;
                imwrite(finalTempMask, [reconstructionDir num2str(nodeItr) '.png']);
            end
        end
        
        % Using the maximum dimensions, transform each composition image to the
        % same size. 
        overallImage = NaN((rowImgCount)*(compMaskSize(1)+1)+1, colImgCount*(compMaskSize(2)+1)+1, size(tempMask,3));
        finalMask = NaN([compMaskSize, size(tempMask,3)]);
        for nodeItr = 1:numberOfNodes
            compFinalMask = finalMask;
            compRealMask = double(allCompMasks{nodeItr});
            margins = (compMaskSize - [size(compRealMask,1), size(compRealMask,2)])/2;
            compFinalMask((floor(margins(1))+1):(end-ceil(margins(1))), ...
                (floor(margins(2))+1):(end-ceil(margins(2))), :) = compRealMask;

            % If this feature is a dead one, reduce illumination by three.
            if isAutoFilter && ismember(nodeItr, deadFeatures)
                compFinalMask = round(compFinalMask * 0.33);
            end
            
            % A make-up to fill in NaNs (empty points).
            fillInValue = median(compFinalMask(~isnan(compFinalMask) & compFinalMask<255));
            compFinalMask(isnan(compFinalMask)) = fillInValue;
            
            % Add the composition's mask to the overall mask image.
            rowStart = 2 + floor((nodeItr-1)/colImgCount)*(compMaskSize(1)+1);
            colStart = 2 + rem(nodeItr-1, colImgCount) * (compMaskSize(2)+1);
            overallImage(rowStart:(rowStart+compMaskSize(1)-1), ...
                colStart:(colStart+compMaskSize(2)-1), :) = compFinalMask;
        end
        
        % A final make up in order to separate masks from each other by 1s.
        overallImage(isnan(overallImage)) = 255;
        overallImage = uint8(overallImage);
        
        % Then, write the compositions the final image.
        imwrite(overallImage, [currentFolder '/debug/' datasetName '/level' num2str(levelId) '_vb.png']);
    end
end



