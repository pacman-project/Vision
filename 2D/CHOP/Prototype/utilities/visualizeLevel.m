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
function [] = visualizeLevel( currentLevel, graphLevel, leafNodes, leafSimilarityMatrix, levelId, ~, numberOfPrevNodes, options, isRedundant)
    % Read options to use in this file.
    currentFolder = options.currentFolder;
    datasetName = options.datasetName;
    useReceptiveField = options.useReceptiveField;
    filtBandCount = size(options.filters{1},3);
    numberOfThreads = options.numberOfThreads;
    childrenPerNode = options.vis.nodeReconstructionChildren;
    if ~useReceptiveField
        % Changed on 18.08.2014 to remove cases where we do not use
        % receptive field. Earlier versions of this file includes code for
        % localization of children in general cases.
        return;
    end
    
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
            if size(mask,1)==1
                mask = imresize(mask, size(mask)*2-3);
            end
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
        patchHalfDims = zeros(numberOfPrevNodes,2);
        lowResponseThrs = zeros(numberOfPrevNodes,1);
        
        for nodeItr = 1:numberOfPrevNodes
            tempImg = double(imread([prevLevelDir num2str(nodeItr) '.png']));
            tempImg = (tempImg - min(min(min(tempImg)))) / (max(max(max(tempImg))) - min(min(min(tempImg))));
            prevNodeMasks(nodeItr) = {tempImg};
            avgPrevNodeMasks(nodeItr) = {mean(tempImg,3)};
            patchHalfDims(nodeItr,:) = [size(prevNodeMasks{nodeItr},1), size(prevNodeMasks{nodeItr},2)]/2;
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
%            parallelVocabNodeSets = cellfun(@(x) currentLevel(x), parallelNodeSets, 'UniformOutput', false);
            parallelVocabNodeSets = parallelNodeSets;
            numberOfThreadsUsed = numberOfThreads;
        else
            parallelNodeSets = {nodeSet};
%            parallelVocabNodeSets = {currentLevel};
            parallelVocabNodeSets = {1:numberOfNodes};
        end
        
        %% Go through each node in the current layer, and reconstuct it to
        % get its mask in the end. Each node is reconstructed using the
        % nodes in the previous layer which contribute to its definition. 
        for setItr = 1:numberOfThreadsUsed
            w = warning('off', 'all');
            nodeSet = parallelNodeSets{setItr};
            vocabNodeSet = parallelVocabNodeSets{setItr};
            
            % Go through each composition in current node set.
            for nodeItr = 1:numel(nodeSet)
                %% Get the children (leaf nodes) from all possible instance in the dataset. Keep the info.
                labelId = vocabNodeSet(nodeItr);
                nodeInstances = find(nodeLabelIds==labelId);
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
                numberOfChildren = numel(children);
                
                %% Give weight the children to get a cleaner representation of the average image for the node.
                % First, we learn base scores, which essentially describe
                % how close each child is to each other. A child does not
                % contribute to itself.
                distMatrix = squareform(pdist(childrenCoords));
                maxDist = max(max(distMatrix));
                scoreMatrix = (maxDist - distMatrix) / maxDist;
                diagIdx = 1:(numberOfChildren+1):numberOfChildren*numberOfChildren;
                scoreMatrix( diagIdx ) = 0;
                
                % Now, weight the scores by the similarity matrix entries.
                weightMatrix = leafSimilarityMatrix(children, children);
                weightMatrix = max(max(weightMatrix)) - weightMatrix;
                
                % Multiply score matrix by weight matrix to get weighted
                % scores (element by element).
                scoreMatrix = scoreMatrix .* weightMatrix;
                scoreMatrix = sum(scoreMatrix,1);
                scoreMatrix = scoreMatrix / max(scoreMatrix);
                
                % Take exponentials of score matrix.
                scoreMatrix = scoreMatrix.^3;
                
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
                    minX = childrenCoords(childItr,1) - patchHalfDims(children(childItr),1);
                    maxX = childrenCoords(childItr,1) + patchHalfDims(children(childItr),1);
                    minY = childrenCoords(childItr,2) - patchHalfDims(children(childItr),2);
                    maxY = childrenCoords(childItr,2) + patchHalfDims(children(childItr),2);
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
                maskMinX = floor(maskMinX);
                maskMinY = floor(maskMinY);
                maskMaxX = ceil(maskMaxX);
                maskMaxY = ceil(maskMaxY);
                childrenCoords = round(childrenCoords - [ones(numel(children),1) * maskMinX, ones(numel(children),1) * maskMinY]);

                %% Write the children's masks to the current mask.
                currentMask = zeros((maskMaxX - maskMinX)+1, (maskMaxY - maskMinY)+1, filtBandCount);
                currentLabelMask = zeros((maskMaxX - maskMinX)+1, (maskMaxY - maskMinY)+1);
                for childItr = 1:numel(children)
                    % Write the child's mask to the output.
                    currentMask((childrenCoords(childItr,1)-floor(patchHalfDims(children(childItr),1))):(childrenCoords(childItr,1)+floor(patchHalfDims(children(childItr),1))), ...
                          (childrenCoords(childItr,2)-floor(patchHalfDims(children(childItr),2))):(childrenCoords(childItr,2)+floor(patchHalfDims(children(childItr),2))),:) = ...
                          prevNodeMasks{children(childItr)}.*scoreMatrix(childItr) + ...
                              currentMask((childrenCoords(childItr,1)-floor(patchHalfDims(children(childItr),1))):(childrenCoords(childItr,1)+floor(patchHalfDims(children(childItr),1))), ...
                                  (childrenCoords(childItr,2)-floor(patchHalfDims(children(childItr),2))):(childrenCoords(childItr,2)+floor(patchHalfDims(children(childItr),2))),:);
                    
                    % Mark label image for the child.
                    currentLabelMask((childrenCoords(childItr,1)-floor(patchHalfDims(children(childItr),1))):(childrenCoords(childItr,1)+floor(patchHalfDims(children(childItr),1))), ...
                         (childrenCoords(childItr,2)-floor(patchHalfDims(children(childItr),2))):(childrenCoords(childItr,2)+floor(patchHalfDims(children(childItr),2)))) = ...
                                 max(currentLabelMask((childrenCoords(childItr,1)-floor(patchHalfDims(children(childItr),1))):(childrenCoords(childItr,1)+floor(patchHalfDims(children(childItr),1))), ...
                         (childrenCoords(childItr,2)-floor(patchHalfDims(children(childItr),2))):(childrenCoords(childItr,2)+floor(patchHalfDims(children(childItr),2)))), ...
                         double(avgPrevNodeMasks{children(childItr)}>lowResponseThrs(children(childItr))) * childItr);
                end
                currentLabelImg = label2rgb(currentLabelMask, 'jet', 'k', 'noshuffle');
                for bandItr = 1:3
                    if filtBandCount>1
                        currentLabelImg(:, :, bandItr) = uint8(double(currentLabelImg(:, :, bandItr)) .* currentMask(:,:,bandItr));
                    else
                        currentLabelImg(:, :, bandItr) = uint8(double(currentLabelImg(:, :, bandItr)) .* currentMask);
                    end
                end
                
                %% If the size of the current mask can be divided by two, pad sides to prevent this.
                dimRems = rem(size(currentMask),2);
                if dimRems(1) == 0
                    currentMask = [currentMask; zeros(1, size(currentMask,2), size(currentMask,3))];
                end
                if dimRems(2) == 0
                    currentMask = [currentMask, zeros(size(currentMask, 1), 1, size(currentMask,3))];
                end
                currentMask = (currentMask - min(min(min(currentMask)))) / (max(max(max(currentMask))) - min(min(min(currentMask))));
                
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
    
    %% Combine all compositions and show them within a single image.
    if ~isRedundant
        % Learn number of rows/columns.
        if levelId == 1
            colImgCount = numberOfNodes;
            rowImgCount = 1;
        else
            colImgCount = ceil(sqrt(numberOfNodes));
            rowImgCount = colImgCount;
        end

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

        % Using the maximum dimensions, transform each composition image to the
        % same size. 
        overallImage = zeros((rowImgCount-1)*compMaskSize(1), colImgCount*compMaskSize(2), size(tempMask,3), 'uint8');
        finalMask = zeros([compMaskSize, size(tempMask,3)], 'uint8');
        for nodeItr = 1:numberOfNodes
            compFinalMask = finalMask;
            compRealMask = allCompMasks{nodeItr};
            margins = (compMaskSize - [size(compRealMask,1), size(compRealMask,2)])/2;
            compFinalMask((floor(margins(1))+1):(end-ceil(margins(1))), ...
                (floor(margins(2))+1):(end-ceil(margins(2))), :) = compRealMask;

            % A small make up. Going to mark sides by adding a line of white
            % padding. Framing each composition.
            compFinalMask([1, end], :, :) = 255;
            compFinalMask(:, [1, end], :) = 255;

            % Add the composition's mask to the overall mask image.
            rowStart = 1 + floor((nodeItr-1)/colImgCount)*compMaskSize(1);
            colStart = 1 + rem(nodeItr-1, colImgCount) * compMaskSize(2);
            overallImage(rowStart:(rowStart+compMaskSize(1)-1), ...
                colStart:(colStart+compMaskSize(2)-1), :) = compFinalMask;
        end

        % Then, write the compositions the final image.
        imwrite(overallImage, [currentFolder '/debug/' datasetName '/level' num2str(levelId) '_vb.png']);
    end
end



