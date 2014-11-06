function [vocabLevel, graphLevel, newDistanceMatrix] = postProcessParts(vocabLevel, graphLevel, nodeDistanceMatrix, options)
    edgeCoords = options.edgeCoords;
    edgeQuantize = options.edgeQuantize/2;
    % Assign new labels of the remaining realizations.
    [remainingComps, ~, IC] = unique([graphLevel.labelId]);
    IC = num2cell(int32(IC));
    [graphLevel.labelId] = deal(IC{:});
    clear IC;
    
    % Eliminate unused compositions from vocabulary.
    vocabLevel = vocabLevel(1, remainingComps);
    newLabelArr = num2cell(int32(1:numel(vocabLevel)));
    [vocabLevel.label] = deal(newLabelArr{:});
    
    %% Find the distance matrix among the remaining parts in vocabLevel.
    edgeCoords((size(edgeCoords,1)+1),:) = [0, 0];
    numberOfNodes = numel(vocabLevel);
    vocabNodeLabels = {vocabLevel.children};
    vocabNodeLabels = cellfun(@(x) double(x), vocabNodeLabels, 'UniformOutput', false);
    vocabEdges = {vocabLevel.adjInfo};
    vocabEdges = cellfun(@(x) double(x), vocabEdges, 'UniformOutput', false);
    newEdge = size(edgeCoords,1);
    vocabNeighborModes = cellfun(@(x) [newEdge; x(:,3)], vocabEdges, 'UniformOutput', false);
    vocabNodePositions = cellfun(@(x) edgeCoords(x,:) - repmat(min(edgeCoords(x,:)), numel(x), 1), vocabNeighborModes, 'UniformOutput', false);
    
    % Sort the nodes inside each vocabulary description.
    vocabSortOrder = cell(size(vocabLevel,1),1);
    for vocabNodeItr = 1:numel(vocabLevel)
        [~, vocabSortOrder{vocabNodeItr}] = sortrows(vocabNodePositions{vocabNodeItr});
    end
    vocabDescriptions = cellfun(@(x,y,z) [x(z)', y(z,:)], vocabNodeLabels, vocabNodePositions, vocabSortOrder, 'UniformOutput', false);
    clear vocabNodeLabels vocabEdges vocabNeighborModes vocabNodePositions vocabSortOrder;
    
    newDistanceMatrix = zeros(numberOfNodes);
    distMatEntries = cell(numberOfNodes,1);
    parfor partItr1 = 1:(numberOfNodes-1)
        description1 = vocabDescriptions{partItr1};
        newEntries = zeros(1, numberOfNodes);
        for partItr2 = (partItr1+1):numberOfNodes; 
            description2 = vocabDescriptions{partItr2}; %#ok<PFBNS>
            matchingCost = InexactMatch(description1, description2, edgeQuantize, nodeDistanceMatrix);
            newEntries(partItr2) = matchingCost;
        end
        distMatEntries(partItr1) = {newEntries};
    end
    if numberOfNodes>1
        newDistanceMatrix(1:(numberOfNodes-1), :) = cat(1, distMatEntries{:});
        newDistanceMatrix = newDistanceMatrix + newDistanceMatrix';
    end
    
    % Normalize distances by the size of compared parts.
    childrenCounts = {vocabLevel.children};
    childrenCounts = cellfun(@(x) numel(x), childrenCounts);
    for partItr = 1:numel(vocabLevel);
        newDistanceMatrix(partItr,:) = newDistanceMatrix(partItr,:) ./ ...
           (max(childrenCounts, repmat(childrenCounts(partItr), 1, numberOfNodes)) * 2 - 1);
    end
    newDistanceMatrix = newDistanceMatrix / max(max(newDistanceMatrix));
end

%> Name: InexactMatch
%>
%> Description: Given two coarse part descriptions, this function tries to
%> match them with the lowest cost possible. Please note that this is
%> essentially graph matching problem, which is NP-Complete. Using very high
%> dimension will result in a drastic performance degradation.
%> Cost of replacing a node: 1
%> Cost of re-positioning a node: distance/maxDistancePossible
%>
%> @param description Description of the first composition. It is of the
%> form: nodeId1 posX posY;
%>       nodeId2 posX posY; ...
%> @param description2 Description of the second composition.
%> @param maxDistance Maximum distance in the node positioning space.
%> @param distanceMatrix NxN matrix with each entry containing the distance
%> between two nodes, indexing that specific entry.
%>
%> @retval lowestCost Minimum matching score.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 06.05.2014
function [lowestCost] = InexactMatch(description, description2, maxDistance, distanceMatrix)
    % Get both descriptions to the same size.
    firstDesSize = size(description,1);
    secDesSize = size(description2,1);
    if firstDesSize > secDesSize
        description2 = [description2; inf(firstDesSize-secDesSize,3)];
    elseif firstDesSize<secDesSize
        description = [description; inf(secDesSize-firstDesSize,3)];
    end
    numberOfChildren = max(firstDesSize, secDesSize);
    
    % Get row permutations of the first description.
    rows = sortrows(perms(1:numberOfChildren));
    
    % Compare each permutation of rows of description to description2. The
    % one which gets the minimum cost is our match.
    lowestCost = inf;
    numberOfRows = size(rows,1);
    for permItr = 1:numberOfRows
        currentCost = 0;
        comparedDescription = description(rows(permItr,:),:);

        % Get valid rows to compare.
        validEdges = ~isinf(comparedDescription(:,1)) & ...
            ~isinf(description2(:,1));
        
        % Estimate node-node distances.
        for nodeItr = 1:numberOfChildren
            if validEdges(nodeItr)
                currentCost = currentCost + distanceMatrix(comparedDescription(nodeItr,1), ...
                                            description2(nodeItr,1));
            end
        end

        % Estimate edge-edge distances.
        currentCost = currentCost + sum(sqrt(sum((comparedDescription(validEdges,2:3) - ...
                     description2(validEdges,2:3)).^2,2)))/maxDistance;
        currentCost = currentCost + numberOfChildren - nnz(validEdges);
        
        % Assign lowest cost if current cost is smaller.
        if currentCost<lowestCost
            lowestCost = currentCost;
            if lowestCost == 0
               break; 
            end
        end
    end
end