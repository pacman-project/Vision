%> Name: visualizeLevel
%>
%> Description: Given the current and previous vocabulary level, this
%> function visualizes the current vocabulary level as a separate patche for
%> each word in the vocabulary level.
%>
%> @param currentLevel Current vocabulary level.
%> @param levelId Identifier of the current level.
%> @param modes Modes of the previous level to reconstruct the features.
%> @param numberOfPrevNodes Number of nodes in previous vocabulary level.
%> @param options Program options.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 10.02.2014
function [] = visualizeLevel( currentLevel, levelId, modes, numberOfPrevNodes, options)
    currentFolder = options.currentFolder;
    datasetName = options.datasetName;
    useReceptiveField = options.useReceptiveField;
    
    numberOfThreads = options.numberOfThreads;
    %% Create output folder structure.
    reconstructionDir = [currentFolder '/debug/' datasetName '/level' num2str(levelId) '/reconstruction/'];
    if ~exist(reconstructionDir, 'dir')
       mkdir(reconstructionDir);
    end
    
    %% In level 1, only print low level filters as first n nodes.
    if levelId == 1
        filterDir = [currentFolder '/filters/' options.filterType '/'];
        numberOfNodes = numel(currentLevel);
        for nodeItr = 1:numberOfNodes
            mask = double(imread([filterDir 'filt' num2str(nodeItr) '.png']));
            mask = imresize(mask, size(mask)*2-3);
            mask = (mask - min(min(mask))) / (max(max(mask)) - min(min(mask)));
            imwrite(mask, [reconstructionDir num2str(nodeItr) '.png']);
        end
    else
        %% In other levels, combine the nodes of previous levels depending on mode info and visualize current level.
        % Read previous layer's masks.
        prevLevelDir = [currentFolder '/debug/' datasetName '/level' num2str(levelId-1) '/reconstruction/'];
        numberOfNodes = numel(currentLevel);
        prevNodeMasks = cell(numberOfPrevNodes,1);
        patchHalfDims = zeros(numberOfPrevNodes,2);
        lowResponseThrs = zeros(numberOfPrevNodes,1);
        
        for nodeItr = 1:numberOfPrevNodes
            tempImg = double(imread([prevLevelDir num2str(nodeItr) '.png']));
            tempImg = (tempImg - min(min(tempImg))) / (max(max(tempImg)) - min(min(tempImg)));
            prevNodeMasks(nodeItr) = {tempImg};
            patchHalfDims(nodeItr,:) = size(prevNodeMasks{nodeItr})/2;
            lowResponseThrs(nodeItr) = max(max(tempImg))/10;
        end
        
        %% To parallelize things, we put vocabulary nodes in different sets, and give each to a thread.
        nodeSet = (1:numberOfNodes);
        numberOfThreadsUsed = 1;
        if options.parallelProcessing && numberOfThreads > 1 && numberOfNodes > numberOfThreads
            parallelNodeSetIdx = round(1:((numberOfNodes-1)/numberOfThreads):numberOfNodes);
            parallelNodeSetIdx = parallelNodeSetIdx(2:end)-parallelNodeSetIdx(1:(end-1));
            parallelNodeSetIdx(end) = parallelNodeSetIdx(end) + 1;
            parallelNodeSets = mat2cell(nodeSet, 1, parallelNodeSetIdx);
            parallelVocabNodeSets = cellfun(@(x) currentLevel(x), parallelNodeSets, 'UniformOutput', false);
            numberOfThreadsUsed = numberOfThreads;
        else
            parallelNodeSets = {nodeSet};
            parallelVocabNodeSets = {currentLevel};
        end
        
        %% Go through each node in the current layer, and reconstuct it to
        % get its mask in the end. Each node is reconstructed using the
        % nodes in the previous layer which contribute to its definition. 
        parfor setItr = 1:numberOfThreadsUsed
            w = warning('off', 'all');
            nodeSet = parallelNodeSets{setItr};
            vocabNodeSet = parallelVocabNodeSets{setItr};
            
            % Go through each composition in current node set.
            for nodeItr = 1:numel(nodeSet)
                node = vocabNodeSet(nodeItr);
                children = node.children;
                childrenCoords = zeros(numel(children),2);

                % Set the first child as placed, such that it is placed into
                % the 0,0 point. All others will be set relatively to that
                % point.
                adjInfo = node.adjInfo;

                %% If receptive fields are used, the children's positions are found with an easier method.
                % If not, a more general localization algorithm is used.
                if useReceptiveField
                    childrenCoords(2:end,:) = -modes(adjInfo(:,3), 3:4);
                else
                    %% Go over the adjacency info array again and again till all nodes are placed.
                    knownNodes = 1;
                    progressing = true;
                    while numel(knownNodes) < numel(children) && progressing
                        progressing = false;
                        edgesToExpand = adjInfo(ismember(adjInfo(:,1), knownNodes) | ...
                            ismember(adjInfo(:,2), knownNodes),:);
                        edgesToExpand = edgesToExpand(~ismember(edgesToExpand(:,1), knownNodes) | ...
                            ~ismember(edgesToExpand(:,2), knownNodes),:);
                        %% Expand each edge. 
                        % This operation basically specifies the points of
                        % previously unknown children. Each edge to expand should
                        % exactly have 1 known node, and 1 unknown node.
                        for edgeItr = 1:size(edgesToExpand,1)
                            if ismember(edgesToExpand(edgeItr,1), knownNodes) && ...
                                ~ismember(edgesToExpand(edgeItr,2), knownNodes)
                                knownNode = edgesToExpand(edgeItr,1);
                                unknownNode = edgesToExpand(edgeItr,2);
                            elseif ismember(edgesToExpand(edgeItr,2), knownNodes) && ...
                                ~ismember(edgesToExpand(edgeItr,1), knownNodes)
                                knownNode = edgesToExpand(edgeItr,2);
                                unknownNode = edgesToExpand(edgeItr,1);
                            else
                                % Do nothing, this node already added in this
                                % iteration.
                                continue;
                            end
                            progressing = true;

                            %% Calculate relative position of the unknown node given known one.
                            modeNo = edgesToExpand(edgeItr,3);
                            isDirected = edgesToExpand(edgeItr,4);
                            relativePosition = modes(modeNo, 3:4);
                            % If node is undirected, determine the positioning by
                            % putting node with smaller type index to the reference
                            % point. On the other hand, if it is directed and known
                            % node is the second one, again take the inverse of
                            % relative position.
                            if (~isDirected && children(knownNode)>children(unknownNode)) || (isDirected && knownNode == edgesToExpand(edgeItr,2))
                                relativePosition = relativePosition * -1;
                            end
                            absolutePosition = relativePosition + childrenCoords(knownNode,:);
                            childrenCoords(unknownNode,:) = absolutePosition;
                            knownNodes = [knownNodes, unknownNode];
                        end
                    end
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
                currentMask = zeros((maskMaxX - maskMinX)+1, (maskMaxY - maskMinY)+1);
                currentLabelMask = zeros((maskMaxX - maskMinX)+1, (maskMaxY - maskMinY)+1);
                for childItr = 1:numel(children)
                      currentMask((childrenCoords(childItr,1)-floor(patchHalfDims(children(childItr),1))):(childrenCoords(childItr,1)+floor(patchHalfDims(children(childItr),1))), ...
                          (childrenCoords(childItr,2)-floor(patchHalfDims(children(childItr),2))):(childrenCoords(childItr,2)+floor(patchHalfDims(children(childItr),2)))) = ...
                          prevNodeMasks{children(childItr)} + ...
                              currentMask((childrenCoords(childItr,1)-floor(patchHalfDims(children(childItr),1))):(childrenCoords(childItr,1)+floor(patchHalfDims(children(childItr),1))), ...
                                  (childrenCoords(childItr,2)-floor(patchHalfDims(children(childItr),2))):(childrenCoords(childItr,2)+floor(patchHalfDims(children(childItr),2))));
%                       currentMask((childrenCoords(childItr,1)-floor(patchHalfDims(children(childItr),1))):(childrenCoords(childItr,1)+floor(patchHalfDims(children(childItr),1))), ...
%                           (childrenCoords(childItr,2)-floor(patchHalfDims(children(childItr),2))):(childrenCoords(childItr,2)+floor(patchHalfDims(children(childItr),2)))) = ...
%                           max(prevNodeMasks{children(childItr)}, ...
%                               currentMask((childrenCoords(childItr,1)-floor(patchHalfDims(children(childItr),1))):(childrenCoords(childItr,1)+floor(patchHalfDims(children(childItr),1))), ...
%                                   (childrenCoords(childItr,2)-floor(patchHalfDims(children(childItr),2))):(childrenCoords(childItr,2)+floor(patchHalfDims(children(childItr),2)))));

    %                             currentMask((childrenCoords(childItr,1)-floor(patchHalfDims(children(childItr),1))):(childrenCoords(childItr,1)+floor(patchHalfDims(children(childItr),1))), ...
    %                     (childrenCoords(childItr,2)-floor(patchHalfDims(children(childItr),2))):(childrenCoords(childItr,2)+floor(patchHalfDims(children(childItr),2)))) = ...
    %                             prevNodeMasks{children(childItr)};
                    currentLabelMask((childrenCoords(childItr,1)-floor(patchHalfDims(children(childItr),1))):(childrenCoords(childItr,1)+floor(patchHalfDims(children(childItr),1))), ...
                         (childrenCoords(childItr,2)-floor(patchHalfDims(children(childItr),2))):(childrenCoords(childItr,2)+floor(patchHalfDims(children(childItr),2)))) = ...
                                 max(currentLabelMask((childrenCoords(childItr,1)-floor(patchHalfDims(children(childItr),1))):(childrenCoords(childItr,1)+floor(patchHalfDims(children(childItr),1))), ...
                         (childrenCoords(childItr,2)-floor(patchHalfDims(children(childItr),2))):(childrenCoords(childItr,2)+floor(patchHalfDims(children(childItr),2)))), ...
                         double(prevNodeMasks{children(childItr)}>lowResponseThrs(children(childItr))) * childItr);
                end
                currentLabelImg = label2rgb(currentLabelMask, 'jet', 'k', 'noshuffle');
                for bandItr = 1:3
                    currentLabelImg(:, :, bandItr) = uint8(double(currentLabelImg(:, :, bandItr)) .* currentMask);
                end
                
                %% If the size of the current mask can be divided by two, pad sides to prevent this.
                dimRems = rem(size(currentMask),2);
                if dimRems(1) == 0
                    currentMask = [currentMask; zeros(1, size(currentMask,2))];
                end
                if dimRems(2) == 0
                    currentMask = [currentMask, zeros(size(currentMask,1), 1)];
                end
                currentMask = (currentMask - min(min(currentMask))) / (max(max(currentMask)) - min(min(currentMask)));
                imwrite(currentMask, [reconstructionDir num2str(nodeSet(nodeItr)) '.png']);
                imwrite(currentLabelImg, [reconstructionDir num2str(nodeSet(nodeItr)) '_comp.png']);
            end
            warning(w);
        end
    end
    
    %% Combine all compositions and show them within a single image.
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



