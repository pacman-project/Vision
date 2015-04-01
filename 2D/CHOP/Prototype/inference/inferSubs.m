%> Name: runSubdue
%>
%> Description: Inference on graphLevel with given subs in vocabLevel. Each
%> composition is searched for in graphLevel.
%>
%> @param vocabLevel Input vocabulary level. Compositions in this vocabulary 
%> level are detected in graphLevel.
%> @param graphLevel The current object graphs' level. The graphLevel's
%> nodes are sorted first by their imageId, then labelId.
%> @param options Program options.
%>
%> @retval graphLevel The graph level consisting of discovered nodes.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 05.02.2014
function [exportArr, confidenceArr] = inferSubs(vocabulary, nodes, distanceMatrices, optimalThresholds, options)
    % Read data into helper data structures.
    scaling = options.scaling;
    edgeQuantize = options.edgeQuantize;
    downEdgeRadius = floor((edgeQuantize-1)/2);
    edgeRadius = options.edgeRadius;
    edgeIdMatrix = options.edgeIdMatrix;
    edgeDistanceMatrix = double(options.edgeDistanceMatrix);
    noveltyThr = 1 - options.noveltyThr;
    halfMatrixSize = (options.edgeQuantize+1)/2;
    matrixSize = [options.edgeQuantize, options.edgeQuantize];
    
    exportArr = [];
    if isempty(nodes) || isempty(vocabulary)
        return;
    end
    nodes = double(nodes);
    allNodes = cell(numel(vocabulary),1);
    allConfidences = cell(numel(vocabulary),1);
    allNodes(1) = {nodes};
    allConfidences(1) = {ones(size(nodes,1), 1, 'double')};
    leafNodeArr = num2cell((1:size(nodes,1))');
    
    for vocabLevelItr = 2:numel(vocabulary)
        scale = (1/options.scaling)^(vocabLevelItr-1);
        neighborhood = floor(options.edgeRadius * scale);
        downsampleRatio = floor((options.edgeQuantize-1)/2) / neighborhood;
 %       scale = (1/options.scaling)^(vocabLevelItr-1);
 %      neighborhood = fix(options.edgeRadius*scale);
        % First, downsample nodes.
%         if vocabLevelItr>2
%             nodes(:,2:3) = round(nodes(:,2:3)) * scaling;
%         end
        
        %% Match subs from vocabLevel to their instance in graphLevel.
        vocabLevel = vocabulary{vocabLevelItr};
        vocabRealizations = cell(numel(vocabLevel),1);
        vocabRealizationsConfidence = cell(numel(vocabLevel),1);
        nodeDistanceMatrix = distanceMatrices{vocabLevelItr-1};
        for vocabItr = 1:numel(vocabLevel)
             adaptiveThreshold = single(optimalThresholds(vocabLevelItr) * ((size(vocabLevel(vocabItr).adjInfo,1)+1)*2-1)) + 0.0001; % Hard threshold for cost of matching two subs.
             %% Subgraph matching.
             % Start with the center.
             centerId = vocabLevel(vocabItr).children(1);
             centerMatchCosts = nodeDistanceMatrix(centerId, nodes(:,1))';
             validInstances = centerMatchCosts < adaptiveThreshold;
             
             if nnz(validInstances) == 0
                 continue;
             end
             numberOfNodes = size(nodes,1);

            % Allocate space to hold the instances.
             instanceChildren = int32(find(validInstances));
             instanceMatchCosts = centerMatchCosts(validInstances);

            %% Get descriptors for edges in the vocabulary node.
            vocabEdges = vocabLevel(vocabItr).adjInfo;

            %% Iteratively, try to match each edge to the instances.
            for edgeItr = 1:size(vocabEdges,1)
                if isempty(instanceChildren) 
                    break;
                end
                childInstances = cell(size(instanceChildren,1),1);
                childMatchingCosts = cell(size(instanceChildren,1),1);
                for parentItr = 1:size(instanceChildren,1)
                    % Find neighbors close to the peripheral node in
                    % question. We're trying to match an edge to the data
                    % in hand.
                    neighborPosDiff = nodes(:,2:3) - repmat(nodes(instanceChildren(parentItr,1),2:3), numberOfNodes, 1);
                    nodeDistances = sqrt(sum(neighborPosDiff.^2, 2));
                    validNeighbors = find(nodeDistances <= neighborhood);
                    validNeighbors = setdiff(validNeighbors, instanceChildren(parentItr,:));
                    
                    % TODO: Check edge novelty, as well as max number of
                   % edges here.
              
                    % Find edge labels.
                    neighborPosDiff = neighborPosDiff(validNeighbors,:);
                    normalizedNeighborPosDiff = fix(fix(downsampleRatio * neighborPosDiff) + halfMatrixSize);
                    matrixInd = sub2ind(matrixSize, normalizedNeighborPosDiff(:,1), normalizedNeighborPosDiff(:,2));
                    neighborEdgeIds = edgeIdMatrix(matrixInd);
                    matchCosts = instanceMatchCosts(parentItr) + edgeDistanceMatrix(neighborEdgeIds, vocabEdges(edgeItr,3)) + ...
                                                                        nodeDistanceMatrix(nodes(validNeighbors,1), vocabLevel(vocabItr).children(vocabEdges(edgeItr,2)));
                                                                    
                    % If any children are found, we save them.
                    validChildrenIdx = matchCosts < adaptiveThreshold;
                    validChildren = validNeighbors(validChildrenIdx);

                    if ~isempty(validChildren)
                        newInstanceChildren = zeros(numel(validChildren), size(instanceChildren,2)+1, 'int32');
                        newInstanceChildren(:,1:size(instanceChildren,2)) = repmat(instanceChildren(parentItr,:), numel(validChildren),1);
                        newInstanceChildren(:,end) = validChildren;
                        childInstances(parentItr) = {newInstanceChildren};
                        childMatchingCosts(parentItr) = {matchCosts(validChildrenIdx)};
                    end
                end
                instanceChildren = cat(1, childInstances{:});
                instanceMatchCosts = cat(1, childMatchingCosts{:});
            end
            
%             numberOfInstances = size(instanceChildren,1);
% %            validInstances = ones(numberOfInstances,1) > 0;
%             % Finally, we eliminate instances that share nodes between them.
%             matchedNodes = zeros(numberOfNodes,1)>0;
%             validInstances = zeros(numberOfInstances,1) > 0;
%             for instItr = 1:numberOfInstances
%                 children = instanceChildren(instItr,:);
%                 if nnz(matchedNodes(instanceChildren(instItr,:))) == 0
%                     matchedNodes(children) = 1;
%                     validInstances(instItr) = 1;
%                 end
%             end

%             instanceChildren = instanceChildren(validInstances,:);
           if isempty(vocabEdges)
               validInstances = zeros(size(instanceChildren,1),1) > 0;
               for parentItr = 1:size(instanceChildren,1)
                    neighborPosDiff = nodes(:,2:3) - repmat(nodes(instanceChildren(parentItr,1),2:3), numberOfNodes, 1);
                    nodeDistances = sum(neighborPosDiff.^2, 2);
                    validNeighbors = find(nodeDistances <= downEdgeRadius^2);
                    validNeighbors = setdiff(validNeighbors, instanceChildren(parentItr,:));
                    if isempty(validNeighbors)
                       validInstances(parentItr) = 1;
                    end
               end
               
               % Eliminate instances which have edges leading out.
               if nnz(validInstances) > 0    
                   instanceChildren = instanceChildren(validInstances, :);
                   instanceMatchCosts = instanceMatchCosts(validInstances, :);
               else
                   instanceChildren = [];
                   instanceMatchCosts = [];
               end
           end

           % Save the instances.
           vocabRealizations(vocabItr) = {instanceChildren};
           instanceConfidences = (adaptiveThreshold - instanceMatchCosts)/ adaptiveThreshold;
           vocabRealizationsConfidence(vocabItr) = {instanceConfidences};
        end

        % Get confidences.
        confidenceArr = double(cat(1, vocabRealizationsConfidence{:}));
        numberOfInstances =numel(confidenceArr);
        
        if numberOfInstances == 0 || vocabLevelItr == numel(vocabulary)
            break;
        end
        
        %% We'll form a new nodes array. 
        % Estimate a new confidence for each realization by propogating
        % from previous layers.
        newNodes = zeros(numberOfInstances,3);
        startIdx = 1;
        for vocabNodeItr = 1:numel(vocabLevel)
           children  = vocabRealizations{vocabNodeItr};
           if isempty(children)
               continue;
           end
           for instanceItr = 1:size(children,1)
               newNodes(startIdx,1) = vocabNodeItr;
               newNodes(startIdx,2:3) = round(mean(nodes(children(instanceItr,:),2:3),1));
               startIdx = startIdx+1;
           end
        end
        
        %% We do two things here in order to eliminate number of nodes.]
        % Step zero: Get the list of leaf nodes for every instance.
        instanceOffset = 1;
        newLeafNodes = cell(numberOfInstances,1);
        for vocabNodeItr = 1:numel(vocabLevel)
            if ~isempty(vocabRealizations{vocabNodeItr})
                children = vocabRealizations{vocabNodeItr};
                for instanceItr = 1:size(children,1)
                    leafNodes = unique(cat(2, leafNodeArr{children(instanceItr,:)}));
                    newLeafNodes{instanceOffset} = leafNodes;
                    instanceOffset = instanceOffset + 1;
                end
            end
        end
        
        % First, we assign each unique set of children to their best matches in the vocabulary.
        [sortedConfidenceArr, sortIdx] = sort(confidenceArr, 'descend');
        sortedLeafNodeArr = newLeafNodes(sortIdx);
        sortedNewNodes = newNodes(sortIdx, :);
        imagePreservedNodes = ones(numberOfInstances,1)>0;
        imageNodeCoords = sortedNewNodes(:,2:3);
        maxSharedLeafNodes = cellfun(@(x) numel(x) * noveltyThr , sortedLeafNodeArr, 'UniformOutput', false);
        for nodeItr = 1:(numberOfInstances-1)
            %% If nobody has erased this node before, it has a right to be in the final graph.
            if imagePreservedNodes(nodeItr) == 0
                continue;
            end
           
            %% Get each neighboring node.
            thisNodeCoords = sortedNewNodes(nodeItr,2:3);
            centerArr = repmat(thisNodeCoords, numberOfInstances, 1);
            distances = sqrt(sum((centerArr - imageNodeCoords).^2, 2));
            adjacentNodes = imagePreservedNodes & distances <= edgeRadius; 
            adjacentNodes(1:nodeItr) = 0;
            selfLeafNodes = sortedLeafNodeArr{nodeItr};

            %% Go over each adjacent node, and apply inhibition if their leaf nodes are too common under current novelty threshold.
            imagePreservedNodes(adjacentNodes) = cellfun(@(x,y) sum(ismembc(x, selfLeafNodes)) <= y, ...
              sortedLeafNodeArr(adjacentNodes), maxSharedLeafNodes(adjacentNodes));
        end
        sortedNewNodes = sortedNewNodes(imagePreservedNodes, :);
        sortedLeafNodeArr = sortedLeafNodeArr(imagePreservedNodes, :);
        sortedConfidenceArr = sortedConfidenceArr(imagePreservedNodes, :);
        
        % Write to output and move on to the next level.
        allNodes(vocabLevelItr) = {sortedNewNodes};
        allConfidences(vocabLevelItr) = {sortedConfidenceArr};
        leafNodeArr = sortedLeafNodeArr;
        nodes = sortedNewNodes;
    end
    confidenceArr = cat(1, allConfidences{:});
    
    numberOfInstances = sum(cellfun(@(x) size(x,1), allNodes));
    %% If no instances have been found, exit.
    if numberOfInstances<1
        exportArr = [];
        return;
    end
    
    %% Save all nodes to an array, and exit.
    exportArr = zeros(numberOfInstances, 5, 'int32');
    exportArr(:,5) = 1;
    startIdx = 1;
    for vocabLevelItr = 1:numel(allNodes)
        thisLevelNodes = allNodes{vocabLevelItr};
        if isempty(thisLevelNodes)
            break;
        end
        exportArr(startIdx:(startIdx + size(thisLevelNodes,1) - 1), 4) = vocabLevelItr;
        exportArr(startIdx:(startIdx + size(thisLevelNodes,1) - 1), 1:3) = thisLevelNodes;
        startIdx = startIdx + size(thisLevelNodes,1);
    end
end

