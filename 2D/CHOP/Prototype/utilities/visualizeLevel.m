%> Name: visualizeLevel
%>
%> Description: Given the current and previous vocabulary level, this
%> function visualizes the current vocabulary level as a separate patche for
%> each word in the vocabulary level.
%>
%> @param currentLevel Current vocabulary level.
%> @param levelId Identifier of the current level.
%> @param modes Modes of the previous level to reconstruct the features.
%> @param currentFolder Current folder.
%> @param options Program options.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 08.12.2013
function [] = visualizeLevel( currentLevel, levelId, modes, currentFolder, options , datasetName)
    %% Create output folder structure.
    reconstructionDir = [currentFolder '/debug/' datasetName '/level' num2str(levelId) '/reconstruction/'];
    if ~exist(reconstructionDir, 'dir')
       mkdir(reconstructionDir);
    end
    
    %% Calculate edge radius.
    scale = (1/options.scaling)^(levelId-1);
    neighborhood = fix(options.edgeRadius * scale);
    
    %% In level 1, only print low level filters as first n nodes.
    if levelId == 1
        filterDir = [currentFolder '/filters/'];
        numberOfNodes = numel(currentLevel);
        for nodeItr = 1:numberOfNodes
            mask = imread([filterDir 'filt' num2str(nodeItr) '.png']);
            imwrite(mask, [reconstructionDir num2str(nodeItr) '.png']);
        end
    else
        %% In other levels, combine the nodes of previous levels depending on mode info and visualize current level.
        % Read previous layer's masks.
        prevLevelDir = [currentFolder '/debug/' datasetName '/level' num2str(levelId-1) '/reconstruction/'];
        numberOfPrevNodes = max(max(modes(:,1:2)));
        numberOfNodes = numel(currentLevel);
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
            node = currentLevel(nodeItr);
            children = node.children;
            childrenCoords = zeros(numel(children),2);
            
            % Set the first child as placed, such that it is placed into
            % the 0,0 point. All others will be set relatively to that
            % point.
            knownNodes = 1;
            adjInfo = node.adjInfo;
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
            currentMask = zeros((maskMaxX - maskMinX)+1, (maskMaxY - maskMinY)+1, 'uint8');
            for childItr = 1:numel(children)
                currentMask((childrenCoords(childItr,1)-floor(patchHalfDims(children(childItr),1))):(childrenCoords(childItr,1)+floor(patchHalfDims(children(childItr),1))), ...
                    (childrenCoords(childItr,2)-floor(patchHalfDims(children(childItr),2))):(childrenCoords(childItr,2)+floor(patchHalfDims(children(childItr),2)))) = ...
                    max(prevNodeMasks{children(childItr)}, ...
                        currentMask((childrenCoords(childItr,1)-floor(patchHalfDims(children(childItr),1))):(childrenCoords(childItr,1)+floor(patchHalfDims(children(childItr),1))), ...
                            (childrenCoords(childItr,2)-floor(patchHalfDims(children(childItr),2))):(childrenCoords(childItr,2)+floor(patchHalfDims(children(childItr),2)))));
            end
            
            %% If the size of the current mask can be divided by two, pad sides to prevent this.
            dimRems = rem(size(currentMask),2);
            if dimRems(1) == 0
                currentMask = [currentMask; zeros(1, size(currentMask,2), 'uint8')];
            end
            if dimRems(2) == 0
                currentMask = [currentMask, zeros(size(currentMask,1), 1, 'uint8')];
            end
            imwrite(currentMask, [reconstructionDir num2str(nodeItr) '.png']);
        end
    end
end

