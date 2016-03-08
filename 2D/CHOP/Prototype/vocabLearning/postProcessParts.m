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
function [vocabLevel, graphLevel, newDistanceMatrix, nodeDistributions] = postProcessParts(vocabLevel, graphLevel, vocabulary, levelItr, options)
    edgeCoords = options.edgeCoords;
    distType = options.distType;
    vocabulary{levelItr} = vocabLevel;
    filterSize = size(options.filters{1});
    halfSize = ceil(filterSize(1)/2); 
    minPixelValue = 1/255;
    if strcmp(options.filterType, 'gabor')
        inhibitionHalfSize = options.gabor.inhibitionRadius;
        stride = options.gabor.stride;
    else
        inhibitionHalfSize = options.auto.inhibitionRadius;
        stride = options.auto.stride;
    end
    halfSearchSize = 2;
    nodeDistributions = [];
    filters = options.filters;
    filters = cellfun(@(x) (x - min(min(x))) / (max(max(x)) - min(min(x))), filters, 'UniformOutput', false);
    
    % As a case study, we replace gabors with 1D gaussian filters
    % stretched.
   filterSize = size(filters{1},1);
    vals = normpdf(1:filterSize, (filterSize+1)/2, 2);
    vals = vals/max(vals);
    firstFilter = repmat(vals, filterSize, 1);
    angle = 180/numel(filters);
    for filterItr = 0:(numel(filters)-1)
         curAngle = -angle * filterItr;
         curFilter = imrotate(firstFilter, curAngle, 'bilinear', 'crop');
         curFilter(curFilter<minPixelValue) = minPixelValue;
         filters{filterItr+1} = curFilter;
    end
    
    visFilters = options.filters;
    visFilters = cellfun(@(x) (x - min(min(x))) / (max(max(x)) - min(min(x))), visFilters, 'UniformOutput', false);
    for filterItr = 1:numel(visFilters)
         curFilter = visFilters{filterItr};
         curFilter(curFilter<minPixelValue) = minPixelValue;
         visFilters{filterItr} = curFilter;
    end
    
    % Assign new labels of the remaining realizations.
    [remainingComps, ~, IC] = unique([graphLevel.labelId]);
    IC = num2cell(int32(IC));
    [graphLevel.labelId] = deal(IC{:});
    [graphLevel.realLabelId] = deal(IC{:});
    clear IC;
    
    % Eliminate unused compositions from vocabulary.
    vocabLevel = vocabLevel(1, remainingComps);
    newLabelArr = num2cell(int32(1:numel(vocabLevel)));
    [vocabLevel.label] = deal(newLabelArr{:});
    vocabulary{levelItr} = vocabLevel;
    
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
   
   % Find the correct image size.
   imageSize = options.receptiveFieldSize * stride * (options.poolDim ^ (levelItr-2)) * (inhibitionHalfSize+1) + halfSize;
   prevImageSize = options.receptiveFieldSize * stride * (options.poolDim ^ (levelItr-3)) * (inhibitionHalfSize+1);
   imageSize = [imageSize, imageSize];
   
   % For now, we make the image bigger.
   imageSize = round(imageSize * 1.5);
   
   %% First, for efficiency, we obtain pixel-level predictions for every part.
   level1Experts = cell(numberOfNodes, 1);
   searchMultiplier = levelItr - 1;
   parfor vocabNodeItr = 1:numberOfNodes
        % Backproject nodes using modal reconstructions.
        nodes = [vocabNodeItr, 0, 0, levelItr];
        experts = projectNode(nodes, vocabulary, 1, 'modal');
%        [~, experts] = optimizeImagination(nodes, vocabulary, imageSize, prevImageSize, filters, visFilters);
     
        % Center the nodes.
        experts = double(experts);
        minX = min(experts(:,2));
        maxX = max(experts(:,2));
        minY = min(experts(:,3));
        maxY = max(experts(:,3));
        midPoint = round([(minX+maxX)/2, (minY+maxY)/2]);
        experts(:,2:3) = experts(:,2:3) - repmat(midPoint, size(experts,1), 1);
        level1Experts{vocabNodeItr} = experts;
   end
   
   % Normalize positions by placing all in the center.
   for vocabNodeItr = 1:numberOfNodes
        experts = level1Experts{vocabNodeItr};
        experts(:,2:3) = round(experts(:,2:3) + repmat(imageSize/2, size(experts,1), 1));
        level1Experts{vocabNodeItr} = experts;
   end
   
   % Comparison of modal reconstructions involves creating a pixel
   % prediction for every pixel, and then looking for matches.
   muImgs = zeros(numberOfNodes, imageSize(1), imageSize(2));
   varImgs = zeros(numberOfNodes, imageSize(1), imageSize(2));
   newDistanceMatrix = zeros(numel(vocabLevel), 'single');

   % Create a blurring filter.
   H = fspecial('gaussian', 5, 1);
   
   % Get product of expert predictions.
   mkdir([options.debugFolder '/level' num2str(levelItr) '/modalProjection/']);
   debugFolder = options.debugFolder;
   parfor vocabNodeItr = 1:numberOfNodes
        [muImg, varImg, ~, ~] = obtainPoE(level1Experts{vocabNodeItr}, [], [], [], imageSize, visFilters);
        muImg = muImg/max(max(muImg));
        blurredMuImg = imfilter(muImg, H, 'replicate');
        muImgs(vocabNodeItr,:,:) = blurredMuImg/max(max(blurredMuImg));
        varImgs(vocabNodeItr,:,:) = varImg;
        imwrite(muImg / max(max(muImg)), [debugFolder '/level' num2str(levelItr) '/modalProjection/' num2str(vocabNodeItr) '.png']);
        imwrite(blurredMuImg/max(max(blurredMuImg)), [debugFolder '/level' num2str(levelItr) '/modalProjection/' num2str(vocabNodeItr) '_blurred.png']);
   end
   
   if strcmp(distType, 'modal')
        % Finally, calculate distances.
        if numel(vocabLevel) > 1
             % Find the distance of two parts using a number of different
             % techniques.
             for vocabNodeItr = 1:(numel(vocabLevel)-1)
                  for vocabNodeItr2 = (vocabNodeItr+1):numel(vocabLevel)
                       distance = findDistance(squeeze(muImgs(vocabNodeItr,:,:)), squeeze(muImgs(vocabNodeItr2,:,:)), ...
                            squeeze(varImgs(vocabNodeItr,:,:)), squeeze(varImgs(vocabNodeItr2,:,:)), distType, halfSearchSize, searchMultiplier);
                       newDistanceMatrix(vocabNodeItr, vocabNodeItr2) = distance;
                       newDistanceMatrix(vocabNodeItr2, vocabNodeItr) = distance;
                  end
             end
        end
   elseif strcmp(distType, 'hog')
        % Histogram of oriented gradients + euclidean distance as distance
        % function.
        sampleDescriptor = HOG(squeeze(muImgs(1,:,:)))';
        descriptors = zeros(numel(vocabLevel), size(sampleDescriptor,2));
        for vocabNodeItr = 1:numel(vocabLevel)
             descriptors(vocabNodeItr,:) = HOG(squeeze(muImgs(vocabNodeItr,:,:)))';
        end
        newDistanceMatrixVect = pdist(descriptors);
        newDistanceMatrix = squareform(newDistanceMatrixVect);
   elseif strcmp(distType, 'hu')
        % Distance by hu moments + euclidean distance.
        descriptors = zeros(numel(vocabLevel), 7);
        for vocabNodeItr = 1:numel(vocabLevel)
             eta_mat = SI_Moment(squeeze(muImgs(vocabNodeItr,:,:))); 
             descriptors(vocabNodeItr,:) = Hu_Moments(eta_mat);
        end
        newDistanceMatrixVect = pdist(descriptors);
        newDistanceMatrix = squareform(newDistanceMatrixVect);
   end
   
   if nnz(newDistanceMatrix) > 0
       newDistanceMatrix = single(newDistanceMatrix/max(max(newDistanceMatrix)));
   end
    %% Finally, we implement the OR nodes here.
    % We apply agglomerative clustering on the generated distance matrix,
    % and group parts based on their similarities. We have a limited number
    % of resources when selecting parts.
    % First, we check for the necessity.
    % All good, group the nodes here.
    newDistanceMatrixVect = squareform(newDistanceMatrix);
    Z = linkage(newDistanceMatrixVect, 'weighted');
    
    %% Find an optimal number of clusters based on Davies-Bouldin index.
    if size(newDistanceMatrix,1) < 3
         clusters = (1:size(newDistanceMatrix,1))';
    else
         clusterStep = 1;
 %        sampleCounts = 2:clusterStep:min([options.reconstruction.numberOfORNodes, ((size(newDistanceMatrix,1))-1), maxIdx]);
         sampleCounts = 2:clusterStep:min([options.reconstruction.numberOfORNodes, ((size(newDistanceMatrix,1))-1)]);
         DBVals = zeros(numel(sampleCounts),1);
         dunnVals = zeros(numel(sampleCounts),1);
         valIndices =  zeros(numel(sampleCounts),1);
         cutoffRatios =  zeros(numel(sampleCounts),1);
         stepItr = 1;
         C = mean(descriptors,1);
         
         for clusterCount  = sampleCounts
 %             maxVal = max(orgMaxVal, Z(end-(clusterCount-2)))-0.00001;
%              clusters = cluster(Z, 'Cutoff', maxVal, 'Criterion', 'distance');
              clusters = cluster(Z, clusterCount);
              
              %% Dunn's index.
              dunnVals(stepItr) = dunns(max(clusters), double(newDistanceMatrix), clusters);
              
              %% Davies-Boulin index.
              % Calculate the index.
              % First, we calculate the cluster centers.
              centers = zeros(clusterCount, size(descriptors,2));
              for centerItr = 1:clusterCount
                   centers(centerItr,:) = mean(descriptors(clusters == centerItr,:),1);
              end
              
              % Within-cluster sum of squares.
              tempMatrix = zeros(size(centers,2), size(centers,2));
              for centerItr = 1:clusterCount
                   idx = clusters == centerItr;
                   vectDiff = descriptors(idx, :) - repmat(centers(centerItr,:), nnz(idx), 1);
                   vals = vectDiff' * vectDiff;
                   tempMatrix = vals + tempMatrix;
              end
              SSW = trace(tempMatrix);
              
              % Between-cluster sum of squares.
              tempMatrix = zeros(clusterCount,size(centers,2));
              for centerItr = 1:clusterCount
                   tempMatrix(centerItr,:) = nnz(clusters == centerItr) * (centers(centerItr,:) - C);
              end
              SSB = tempMatrix' * tempMatrix;
              SSB = trace(SSB);
              
              % Total sum of squares.
%               vectDiff = descriptors - repmat(C, size(descriptors,1), 1);
%               SST = vectDiff' * vectDiff;
%               SST = trace(SST);
              SST = SSW + SSB;
              
              % Intra-cluster distance.
              tempMatrix = zeros(clusterCount,1);
              for centerItr = 1:clusterCount
                   idx = clusters == centerItr;
                   vectDiff = descriptors(idx, :) - repmat(centers(centerItr,:), nnz(idx), 1);
                   sums = sum(sum(vectDiff.^2));
                   tempMatrix(centerItr) = sqrt(sums);
              end
              CIntra = sum(tempMatrix);
              
              % Inter-cluster distance.
              tempMatrix = zeros(clusterCount, clusterCount);
              for centerItr1 = 1:clusterCount
                   for centerItr2 = 1:clusterCount
                        tempMatrix(centerItr1, centerItr2) = sqrt(sum((centers(centerItr1, :) - centers(centerItr2,:)).^2));
                   end
              end
              CInter = sum(sum(tempMatrix)) / (clusterCount^2);
%              CInter2 = sum(sum(pdist(centers))) / (clusterCount^2);

              % Final index
              valIndices(stepItr) = abs(((SSW/SSB)*SST) - (CIntra/CInter) - (size(descriptors,1) - clusterCount));
              if stepItr > 1
                   if valIndices(stepItr)>valIndices(stepItr-1)
                        cutoffRatios(stepItr) = valIndices(stepItr-1) / valIndices(stepItr);
                   else
                        cutoffRatios(stepItr) = valIndices(stepItr) / valIndices(stepItr-1);
                   end
              end
              
              % Second, We obtain the average distance of cluster samples
              % within every cluster to the center.
              sigmas = zeros(clusterCount, 1);
              for centerItr = 1:clusterCount
                   idx = clusters == centerItr;
                   vectDiff = descriptors(idx, :) - repmat(centers(centerItr,:), nnz(idx), 1);
                   sigmas(centerItr) = mean(sqrt(sum(vectDiff.^2,2)));
              end
              
              % Finally, calculate DB index.
              DB = 0;
              for centerItr = 1:clusterCount
                   tempArr = zeros(clusterCount,1);
                   for centerItr2 = 1:clusterCount
                        if centerItr == centerItr2
                             continue;
                        end
                        tempArr(centerItr2) = (sigmas(centerItr) + sigmas(centerItr2)) / sqrt(sum((centers(centerItr,:) - centers(centerItr2,:)).^2));
                   end
                   DB = DB + max(tempArr);
              end
              DB = DB / clusterCount;
              DBVals(stepItr) = DB;
              stepItr = stepItr + 1;
         end
         
         % Show DB index.
         figure, plot(sampleCounts, cutoffRatios);
         saveas(gcf, [options.debugFolder '/level' num2str(levelItr) '_ValidityIdx.png']);
         close all;
         figure, plot(sampleCounts, DBVals);
         saveas(gcf, [options.debugFolder '/level' num2str(levelItr) '_DBIdx_NoNorm.png']);
         close all;
         figure, plot(sampleCounts, dunnVals);
         saveas(gcf, [options.debugFolder '/level' num2str(levelItr) '_DunnIndex.png']);
         stepSize = 5;
         if size(newDistanceMatrix,1) > 1
              ZVals = Z(2:end,3) - Z(1:(end-1),3);
              figure, plot(1:size(Z,1), Z(:,3));
              saveas(gcf, [options.debugFolder '/level' num2str(levelItr) '_combinationCosts.png']);         
              close all;
              figure, plot(1:numel(ZVals), ZVals);
              saveas(gcf, [options.debugFolder '/level' num2str(levelItr) '_addedCosts.png']);
              close all;
              idx = stepSize:stepSize:size(Z,1);
              figure, plot(idx, Z(idx));
              saveas(gcf, [options.debugFolder '/level' num2str(levelItr) '_combinationCosts_Steps.png']);
              close all;
         end
         
         % Find an ideal cutoff ratio.
         bufSize = 10;
         val = find(cutoffRatios(bufSize:end-bufSize) <= 0.5, 1, 'first');
         clusterCount = sampleCounts(val);
         if isempty(val)
              [~, idx] = max(dunnVals);
              clusterCount = idx(1);
         end
         
         % Find optimal number of clusters.
%          [~, maxIdx] = max(dunnVals);
%          clusterCount = sampleCounts(maxIdx);
         clusters = cluster(Z, 'maxclust', clusterCount);

         % Visualize dendogram.
          figure, dendrogram(Z, min(clusterCount, 50));
          saveas(gcf, [options.debugFolder '/level' num2str(levelItr) '_dendogram.png']);
          close all;
    end
    
    
    % Combine parts falling in the same clusters by setting their distances to zero.
    [~, IA, IC] = unique(clusters, 'stable');
    clusterCount = numel(IA);
    condensedDistanceMatrix = zeros(clusterCount, 'single');
    if numel(unique(clusters)) ~= size(newDistanceMatrix,1)
         for clusterItr1 = 1:clusterCount
              for clusterItr2 = (clusterItr1+1):clusterCount
                   cluster1Samples = IC==clusterItr1;
                   cluster2Samples = IC==clusterItr2;
                   distanceMatrixSmall = newDistanceMatrix(cluster1Samples, cluster2Samples);
                   avgDistance = mean(distanceMatrixSmall(:));
                   condensedDistanceMatrix(clusterItr1, clusterItr2) = avgDistance;
                   condensedDistanceMatrix(clusterItr2, clusterItr1) = avgDistance;
              end
         end
         
         labelArr = num2cell(int32(IC));
         % Update the labels of vocabLevel with or node information.
         [vocabLevel.label] = deal(labelArr{:});
         newDistanceMatrix = condensedDistanceMatrix;
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