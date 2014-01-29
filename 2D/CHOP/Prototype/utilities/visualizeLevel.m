%> Name: visualizeLevel
%>
%> Description: Given the current and previous vocabulary level, this
%> function visualizes the current vocabulary level as a separate patche for
%> each word in the vocabulary level.
%>
%> @param currentVocabLevel Current vocabulary level.
%> @param currentGraphLevel Current graph level.
%> @param levelId Identifier of the current level.
%> @param modes Modes of the previous level to reconstruct the features.
%> @param options Program options.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 08.12.2013
function [] = visualizeLevel( currentVocabLevel, mainGraph, levelId, options)
    %% Create output folder structure.
    reconstructionDir = [options.currentFolder '/debug/' options.datasetName '/level' num2str(levelId) '/reconstruction/'];
    if ~exist(reconstructionDir, 'dir')
       mkdir(reconstructionDir);
    end
    
    %% In level 1, only print low level filters as first n nodes.
    if levelId == 1
        filterDir = [options.currentFolder '/filters/' options.filterType '/'];
        numberOfNodes = numel(currentVocabLevel);
        for nodeItr = 1:numberOfNodes
            mask = imread([filterDir 'filt' num2str(nodeItr) '.png']);
            imwrite(mask, [reconstructionDir num2str(nodeItr) '.png']);
        end
    else
        currentGraphLevel = mainGraph{levelId};
        prevGraphLevel = mainGraph{levelId-1};
        prevGraphLevelIds = [prevGraphLevel.labelId];
        %% In other levels, combine the nodes of previous levels and visualize current level.
        % Read previous layer's masks.
        prevLevelDir = [options.currentFolder '/debug/' options.datasetName '/level' num2str(levelId-1) '/reconstruction/'];
        numberOfPrevNodes = max(prevGraphLevelIds);
        numberOfNodes = numel(currentVocabLevel);
        prevNodeMasks = cell(numberOfPrevNodes,1);
        patchHalfDims = zeros(numberOfPrevNodes,2);
        for nodeItr = 1:numberOfPrevNodes
            prevNodeMasks(nodeItr) = {imread([prevLevelDir num2str(nodeItr) '.png'])};
            patchHalfDims(nodeItr,:) = size(prevNodeMasks{nodeItr})/2;
        end
        
        % Go through each node in the current layer, and reconstuct it to
        % get its mask in the end. Each node is reconstructed using the
        % nodes in the previous layer which contribute to its definition. 
        for nodeItr = 1:numberOfNodes
            %% Access all realizations of this composition and write each to a
            % different mask. The resulting mask corresponding to the
            % composition is basically their average image.
            realizations = currentGraphLevel(1,[currentGraphLevel.labelId]==nodeItr);
            maxMaskSize = [0,0];
            realMasks = cell(numel(realizations),1);
            %% Go over the realizations and write each to a separate mask.
            for realItr = 1:numel(realizations)
                node = realizations(realItr);
                children = node.children;
                childrenCoords = zeros(numel(children),2);

                % Set the first child as placed, such that it is placed into
                % the 0,0 point. All others will be set relatively to that
                % point.
                knownNodes = children(1);
                adjInfo = node.childrenAdjInfo;
                progressing = true;

                % Get only half of directed edges, and all undirected
                % edges.
    %            sufficientAdjInfo = adjInfo((adjInfo(:,1)<adjInfo(:,2) & adjInfo(:,4)==1) | adjInfo(:,4)==0, :);

                %% Go over the adjacency info array again and again till all nodes are placed.
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
        %                modeNo = edgesToExpand(edgeItr,3);
        %                isDirected = edgesToExpand(edgeItr,4);
        %                relativePosition = modes(modeNo, 3:4);
                        position1 = prevGraphLevel(knownNode).position;
                        position2 = prevGraphLevel(unknownNode).position;
                        relativePosition = position2-position1;
        
                        % If node is undirected, determine the positioning by
                        % putting node with smaller type index to the reference
                        % point. On the other hand, if it is directed and known
                        % node is the second one, again take the inverse of
                        % relative position.
 %                       if (~isDirected && prevGraphLevelIds(children(knownNode))>prevGraphLevelIds(children(unknownNode))) || ...
 %                               (isDirected && knownNode == edgesToExpand(edgeItr,2))
 %                           relativePosition = relativePosition * -1;
 %                       end
                        knownNodeIdx = find(children==knownNode, 1, 'first');
                        unknownNodeIdx = find(children==unknownNode, 1, 'first');
                        absolutePosition = relativePosition + childrenCoords(knownNodeIdx,:);
                        childrenCoords(unknownNodeIdx,:) = absolutePosition;
                        knownNodes = [knownNodes, unknownNode];
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
                children = prevGraphLevelIds(children);
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

                %% ADD the children's mask to the current mask.
                currentMask = zeros((maskMaxX - maskMinX)+1, (maskMaxY - maskMinY)+1, 'uint8');
                for childItr = 1:numel(children)
                     currentMask((childrenCoords(childItr,1)-floor(patchHalfDims(children(childItr),1))):(childrenCoords(childItr,1)+floor(patchHalfDims(children(childItr),1))), ...
                         (childrenCoords(childItr,2)-floor(patchHalfDims(children(childItr),2))):(childrenCoords(childItr,2)+floor(patchHalfDims(children(childItr),2)))) = ...
                         max(prevNodeMasks{children(childItr)}, ...
                             currentMask((childrenCoords(childItr,1)-floor(patchHalfDims(children(childItr),1))):(childrenCoords(childItr,1)+floor(patchHalfDims(children(childItr),1))), ...
                                 (childrenCoords(childItr,2)-floor(patchHalfDims(children(childItr),2))):(childrenCoords(childItr,2)+floor(patchHalfDims(children(childItr),2)))));
    %                            currentMask((childrenCoords(childItr,1)-floor(patchHalfDims(children(childItr),1))):(childrenCoords(childItr,1)+floor(patchHalfDims(children(childItr),1))), ...
    %                    (childrenCoords(childItr,2)-floor(patchHalfDims(children(childItr),2))):(childrenCoords(childItr,2)+floor(patchHalfDims(children(childItr),2)))) = ...
    %                            prevNodeMasks{children(childItr)};
                end

                %% If the size of the current mask can be divided by two, pad sides to prevent this.
                dimRems = rem(size(currentMask),2);
                if dimRems(1) == 0
                    currentMask = currentMask(1:(end-1),:);
                end
                if dimRems(2) == 0
                    currentMask = currentMask(:,1:(end-1));
                end
                realMasks{realItr} = currentMask;
                
                % Save mask size if exceeds current boundaries.
                if maxMaskSize(1) < size(currentMask,1)
                   maxMaskSize(1) = size(currentMask,1); 
                end
                if maxMaskSize(2) < size(currentMask,2)
                   maxMaskSize(2) = size(currentMask,2); 
                end
                
         %       [reconstructionDir num2str(nodeItr) '.png'];
            end
            %% Finally, average all resulting images.
            finalCompMask = zeros(maxMaskSize);
            for realItr = 1:numel(realizations)
                patchSize = size(realMasks{realItr});
                padding = ceil((maxMaskSize-patchSize)/2);
                finalCompMask((1+padding(1)):(end-padding(1)), (1+padding(2)):(end-padding(2)), :) = ...
                   finalCompMask((1+padding(1)):(end-padding(1)), (1+padding(2)):(end-padding(2)), :) + ...
                        double(realMasks{realItr});
            end
            finalCompMask = finalCompMask-min(min(finalCompMask));
            finalCompMask = uint8((finalCompMask/max(max(finalCompMask))) * 255);
            imwrite(finalCompMask, [reconstructionDir num2str(nodeItr) '.png']);
        end
    end
end

