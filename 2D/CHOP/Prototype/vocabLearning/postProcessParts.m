%> Name: postProcessParts
%>
%> Description: Finds the distance between the parts in vocabLevel.
%> Additionally, it assumes that we have already removed/inhibited some nodes
%> from graphLevel (object graphs), so it removes parts which do not have any
%> instances from vocabLevel. Once they are removed, the parts are ordered by
%> the number of occurences, while the original ordering is preserved in
%> graphLabelAssgnArr to be used in inference. 
%>
%> @param vocabLevel Vocabulary level to be processed.
%> @param graphLevel Object graphs, encoded in a node + adjacency list
%> fashion.
%> @param nodeDistanceMatrix The distance matrix of previous layer's nodes.
%> @param options Program options.
%>
%> @retval vocabLevel Processed vocabulary level.
%> @retval graphLevel Remaining nodes for object graphs, encoded in a 
%> node + adjacency list fashion.
%> @retval newDistanceMatrix Generated distance matrix of this layer's nodes.
%> @retval graphLabelAssgnArr Original ordering of parts in vocabLevel 
%> before occurence-based reordering.
%>
%> @retval lowestCost Minimum matching score.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 06.05.2014
%> Update on 23.02.2015 Added comments, performance boost.
%> Update on 25.02.2015 Added support for single node subs.
function [vocabLevel, graphLevel, newDistanceMatrix] = postProcessParts(vocabLevel, graphLevel, vocabulary, levelItr, options)
    edgeCoords = options.edgeCoords;
    distType = options.distType;
    vocabulary{levelItr} = vocabLevel;
    filterSize = size(options.filters{1});
    halfSize = ceil(filterSize(1)/2);
    halfSearchSize = 1;
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
   vocabEdges = {vocabLevel.adjInfo};
   vocabEdges = cellfun(@(x) double(x), vocabEdges, 'UniformOutput', false);
   newEdge = size(edgeCoords,1);
   largeSubIdx = cellfun(@(x) ~isempty(x), vocabEdges);
   vocabNeighborModes = num2cell(repmat(newEdge, 1, numel(vocabEdges)));
   vocabNeighborModes(largeSubIdx) = cellfun(@(x,y) [y; x(:,3)], vocabEdges(largeSubIdx), vocabNeighborModes(largeSubIdx), 'UniformOutput', false);
   vocabNodePositions = cellfun(@(x) edgeCoords(x,:) - repmat(min(edgeCoords(x,:)), numel(x), 1), vocabNeighborModes, 'UniformOutput', false);

   % Sort the nodes inside each vocabulary description.
   vocabSortOrder = cell(size(vocabLevel,1),1);
   for vocabNodeItr = 1:numel(vocabLevel)
       [~, vocabSortOrder{vocabNodeItr}] = sortrows(vocabNodePositions{vocabNodeItr});
   end
   clear vocabNodeLabels vocabEdges vocabNeighborModes vocabNodePositions vocabSortOrder;

   %% We  are experimenting with different distance functions.
   if strcmp(distType, 'modal')
        %% First, for efficiency, we obtain pixel-level predictions for every part.
        level1Experts = cell(numberOfNodes, 1);
        minX = 0;
        maxX = 0;
        minY = 0;
        maxY = 0;
        for vocabNodeItr = 1:numberOfNodes
             % Backproject nodes using modal reconstructions.
             nodes = [vocabNodeItr, 0, 0, levelItr];
             experts = projectNode(nodes, vocabulary, 1, 'modal');
             
             % Center the nodes.
             experts = double(experts);
             experts(:,2:3) = experts(:,2:3) - repmat(round(mean(experts(:,2:3),1)), size(experts,1), 1);
             
             % Obtain min/max x,y values.
             minX = min(minX, min(experts(:,2)) - (halfSize + halfSearchSize));
             maxX = max(maxX, max(experts(:,2)) + (halfSize + halfSearchSize));
             minY = min(minY, min(experts(:,3)) - (halfSize + halfSearchSize));
             maxY = max(maxY, max(experts(:,3)) + (halfSize + halfSearchSize));
             level1Experts{vocabNodeItr} = experts;
        end
        
        % Find the correct image size.
        imageSize = [maxX - minX + 1, maxY - minY + 1];
        
        % Normalize positions by placing all in the center.
        for vocabNodeItr = 1:numberOfNodes
             experts = level1Experts{vocabNodeItr};
             experts(:,2:3) = experts(:,2:3) + repmat(ceil(imageSize/2), size(experts,1), 1);
             level1Experts{vocabNodeItr} = experts;
        end
        
        % Comparison of modal reconstructions involves creating a pixel
        % prediction for every pixel, and then looking for matches.
        muImgs = zeros(numberOfNodes, imageSize(1), imageSize(2));
        varImgs = zeros(numberOfNodes, imageSize(1), imageSize(2));
        newDistanceMatrix = zeros(numel(vocabLevel), 'single');
        
        % Get product of expert predictions.
        for vocabNodeItr = 1:numberOfNodes
             [muImg, varImg] = obtainPoE(level1Experts{vocabNodeItr}, imageSize, options);
             muImgs(vocabNodeItr,:,:) = muImg;
             varImgs(vocabNodeItr,:,:) = varImg;
        end
        
        % Finally, calculate distances.
        if numel(vocabLevel) > 1
             % Find the distance of two parts using a number of different
             % techniques.
             for vocabNodeItr = 1:(numel(vocabLevel)-1)
                  for vocabNodeItr2 = (vocabNodeItr+1):numel(vocabLevel)
                       distance = findDistance(squeeze(muImgs(vocabNodeItr,:,:)), squeeze(muImgs(vocabNodeItr2,:,:)), ...
                            squeeze(varImgs(vocabNodeItr,:,:)), squeeze(varImgs(vocabNodeItr2,:,:)), distType, halfSearchSize);
                       newDistanceMatrix(vocabNodeItr, vocabNodeItr2) = distance;
                       newDistanceMatrix(vocabNodeItr2, vocabNodeItr) = distance;
                  end
             end
        end
   end
   
   
   
    %% Finally, we implement the OR nodes here.
    % We apply agglomerative clustering on the generated distance matrix,
    % and group parts based on their similarities. We have a limited number
    % of resources when selecting parts.
    % First, we check for the necessity.
    
    % All good, group the nodes here.
    newDistanceMatrixVect = squareform(newDistanceMatrix);
    Z = linkage(newDistanceMatrixVect, 'average');
    clusters = cluster(Z, 'maxclust', options.reconstruction.numberOfORNodes);
    
    % Combine parts falling in the same clusters by setting their distances to zero.
    [~, IA, IC] = unique(clusters, 'stable');
    numberOfClusters = numel(IA);
    condensedDistanceMatrix = zeros(numberOfClusters, 'single');
    if numel(unique(clusters)) ~= size(newDistanceMatrix,1)
         for clusterItr1 = 1:numberOfClusters
              for clusterItr2 = (clusterItr1+1):numberOfClusters
                   cluster1Samples = IC==clusterItr1;
                   cluster2Samples = IC==clusterItr2;
                   distanceMatrixSmall = newDistanceMatrix(cluster1Samples, cluster2Samples);
                   avgDistance = mean(distanceMatrixSmall(:));
                   condensedDistanceMatrix(clusterItr1, clusterItr2) = avgDistance;
                   condensedDistanceMatrix(clusterItr2, clusterItr1) = avgDistance;
              end
         end
    end
    labelArr = num2cell(int32(IC));
    
    % Update the labels of vocabLevel with or node information.
    [vocabLevel.label] = deal(labelArr{:});
    newDistanceMatrix = condensedDistanceMatrix;
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
    lowestCost = realmax('single');
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
                currentCost = currentCost + single(full(distanceMatrix(comparedDescription(nodeItr,1), ...
                                            description2(nodeItr,1))));
            else
                 currentCost = currentCost + 1;
            end
        end

        % Estimate edge-edge distances.
        currentCost = currentCost + sum(sqrt(sum((comparedDescription(validEdges,2:3) - ...
                     description2(validEdges,2:3)).^2,2)))/maxDistance;
        currentCost = currentCost + numberOfChildren - nnz(validEdges);
        
        % Assign lowest cost if current cost is smaller.
        if currentCost<lowestCost
            lowestCost = currentCost;
        end
    end
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
function [lowestCost] = CalculateLogLikelihood(description, description2, maxDistance, distanceMatrix)
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
    lowestCost = realmax('single');
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
                currentCost = currentCost + single(full(distanceMatrix(comparedDescription(nodeItr,1), ...
                                            description2(nodeItr,1))));
            else
                 currentCost = currentCost + 1;
            end
        end

        % Estimate edge-edge distances.
        currentCost = currentCost + sum(sqrt(sum((comparedDescription(validEdges,2:3) - ...
                     description2(validEdges,2:3)).^2,2)))/maxDistance;
        currentCost = currentCost + numberOfChildren - nnz(validEdges);
        
        % Assign lowest cost if current cost is smaller.
        if currentCost<lowestCost
            lowestCost = currentCost;
        end
    end
end