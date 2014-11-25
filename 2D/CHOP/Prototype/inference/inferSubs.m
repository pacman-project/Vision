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
function [exportArr] = inferSubs(vocabulary, nodes, distanceMatrices, options)
    % Read data into helper data structures.
    threshold = options.subdue.threshold;
    scaling = options.scaling;
    edgeQuantize = options.edgeQuantize;
    downEdgeRadius = floor((edgeQuantize-1)/2);
    edgeRadius = options.edgeRadius;
    edgeIdMatrix = options.edgeIdMatrix;
    edgeDistanceMatrix = double(options.edgeDistanceMatrix);
    
    exportArr = [];
    if isempty(nodes) || isempty(vocabulary)
        return;
    end
    nodes = double(nodes);
    nodes(:,2:3) = (nodes(:,2:3)) * (downEdgeRadius / edgeRadius);
    allNodes = cell(numel(vocabulary),1);
    allNodes(1) = {nodes};
    for vocabLevelItr = 2:numel(vocabulary)
        % First, downsample nodes.
        if vocabLevelItr>2
            nodes(:,2:3) = round(nodes(:,2:3)) * scaling;
        end
        
        %% Match subs from vocabLevel to their instance in graphLevel.
        vocabLevel = vocabulary{vocabLevelItr};
        vocabRealizations = cell(numel(vocabLevel),1);
        vocabRealizationsConfidence = cell(numel(vocabLevel),1);
        nodeDistanceMatrix = distanceMatrices{vocabLevelItr-1};
        for vocabItr = 1:numel(vocabLevel)
             adaptiveThreshold = single(threshold * ((size(vocabLevel(vocabItr).adjInfo,1)+1)*2-1)) + 0.0001; % Hard threshold for cost of matching two subs.
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
                    nodeDistances = sum(neighborPosDiff.^2, 2);
                    validNeighbors = find(nodeDistances <= downEdgeRadius^2);
                    validNeighbors = setdiff(validNeighbors, instanceChildren(parentItr,1));
                    neighborPosDiff = neighborPosDiff(validNeighbors,:);
                    edgeIdPos = fix(neighborPosDiff) + downEdgeRadius + 1;
                    neighborEdgeIdx = sub2ind([edgeQuantize, edgeQuantize], edgeIdPos(:,1), edgeIdPos(:,2));
                    neighborEdgeIds = edgeIdMatrix(neighborEdgeIdx);
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
           vocabRealizations(vocabItr) = {instanceChildren};
           instanceConfidences = (adaptiveThreshold - instanceMatchCosts )/ adaptiveThreshold;
           vocabRealizationsConfidence(vocabItr) = {instanceConfidences};
        end

        % Get confidences.
        confidenceArr = double(cat(1, vocabRealizationsConfidence{:}));
        confidenceArrCell = num2cell(confidenceArr);
        numberOfInstances =numel(confidenceArr);
        
        if numberOfInstances == 0 || vocabLevelItr == numel(vocabulary)
            break;
        end
        
        %% Here, we'll form a new nodes array. 
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
        allNodes(vocabLevelItr) = {newNodes};
        nodes = newNodes;
    end
    
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
        allNodes(startIdx:(startIdx + numel(thisLevelNodes) - 1), 4) = vocabLevelItr;
        allNodes(startIdx:(startIdx + numel(thisLevelNodes) - 1), 1:3) = thisLevelNodes;
    end
end

