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
    stopFlag = options.stopFlag;
    
    vocabulary{levelItr} = vocabLevel;
    filterSize = size(options.filters{1});
    halfSize = ceil(filterSize(1)/2); 
    if strcmp(options.filterType, 'gabor')
        inhibitionHalfSize = options.gabor.inhibitionRadius;
        stride = options.gabor.stride;
    else
        inhibitionHalfSize = options.auto.inhibitionRadius;
        stride = options.auto.stride;
    end
    
    %% As a case study, we replace gabors with 1D gaussian filters    
    visFilters = options.filters;
    visFilters = cellfun(@(x) (x - min(min(x))) / (max(max(x)) - min(min(x))), visFilters, 'UniformOutput', false);
    for filterItr = 1:numel(visFilters)
         visFilters{filterItr} = uint8(round(visFilters{filterItr} * 255));
    end
    visFilters = cat(3, visFilters{:});
    visFilters(visFilters<1) = 1;
    
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
   imageSize = [imageSize, imageSize];
   
   % For now, we make the image bigger.
%   imageSize = round(imageSize * 1.5);
   
   %% First, for efficiency, we obtain pixel-level predictions for every part.
   level1Experts = cell(numberOfNodes, 1);
   parfor vocabNodeItr = 1:numberOfNodes
        % Backproject nodes using modal reconstructions.
        nodes = [vocabNodeItr, 0, 0, levelItr];
        experts = projectNode(nodes, vocabulary, 1, 'modal');

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
   newDistanceMatrix = zeros(numel(vocabLevel), 'single');

   % Create a blurring filter.
   H = fspecial('gaussian', 5, 1);
   
   % Get product of expert predictions.
   display('........ Visualizing parts using Product of Experts..');
   mkdir([options.debugFolder '/level' num2str(levelItr) '/modalProjection/']);
   debugFolder = options.debugFolder;
   parfor vocabNodeItr = 1:numberOfNodes
        [muImg, ~, ~] = obtainPoE(level1Experts{vocabNodeItr}, [], [], imageSize, visFilters, []);
        muImg = uint8(round(255*(double(muImg) / double(max(max(muImg))))));
        blurredMuImg = uint8(imfilter(muImg, H, 'replicate'));
        muImgs(vocabNodeItr,:,:) = blurredMuImg;
        blurredMuImg = uint8(round(255*(double(blurredMuImg) / double(max(max(blurredMuImg))))));
        imwrite(muImg, [debugFolder '/level' num2str(levelItr) '/modalProjection/' num2str(vocabNodeItr) '.png']);
        imwrite(blurredMuImg, [debugFolder '/level' num2str(levelItr) '/modalProjection/' num2str(vocabNodeItr) '_blurred.png']);
   end
   
   display('........ Calculating part distances and performing clustering, based on a dynamic cut-off index..');
   if strcmp(distType, 'modal')
        % Finally, calculate distances.
        % TODO: Currently doesn't work! We have removed variance output, we
        % need to remove that from findDistance as well.
%         if numel(vocabLevel) > 1
%              % Find the distance of two parts using a number of different
%              % techniques.
%              for vocabNodeItr = 1:(numel(vocabLevel)-1)
%                   for vocabNodeItr2 = (vocabNodeItr+1):numel(vocabLevel)
%                        distance = findDistance(squeeze(muImgs(vocabNodeItr,:,:)), squeeze(muImgs(vocabNodeItr2,:,:)), ...
%                             squeeze(varImgs(vocabNodeItr,:,:)), squeeze(varImgs(vocabNodeItr2,:,:)), distType, halfSearchSize, searchMultiplier);
%                        newDistanceMatrix(vocabNodeItr, vocabNodeItr2) = distance;
%                        newDistanceMatrix(vocabNodeItr2, vocabNodeItr) = distance;
%                   end
%              end
%         end
   elseif strcmp(distType, 'hog')
        % Histogram of oriented gradients + euclidean distance as distance
        % function.
        sampleDescriptor = HOG(squeeze(muImgs(1,:,:)))';
        descriptors = zeros(numel(vocabLevel), size(sampleDescriptor,2));
        parfor vocabNodeItr = 1:numel(vocabLevel)
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
    Z = linkage(newDistanceMatrixVect, 'ward');
    
    %% Find an optimal number of clusters based on Davies-Bouldin index.
    if size(newDistanceMatrix,1) < 3
         clusters = (1:size(newDistanceMatrix,1))';
    else
         % If stop flag is not up, we find an optimal number of clusters.
         if ~stopFlag
              clusterStep = 1;
              maxNumberOfORNodes = options.reconstruction.maxNumberOfORNodes;
      %        sampleCounts = 2:clusterStep:min([options.reconstruction.numberOfORNodes, ((size(newDistanceMatrix,1))-1), maxIdx]);
              sampleCounts = min([round(maxNumberOfORNodes/10), ((size(newDistanceMatrix,1))-1)]):clusterStep:min([maxNumberOfORNodes, ((size(newDistanceMatrix,1))-1)]);
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
     %              SSWTemp = 0;
                   for centerItr = 1:clusterCount
                        idx = clusters == centerItr;
                        vectDiff = descriptors(idx, :) - repmat(centers(centerItr,:), nnz(idx), 1);
     %                   SSWCluster = sum(sum(vectDiff.^2));
     %                   SSWTemp = SSWTemp + SSWCluster;
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
                    vectDiff = descriptors - repmat(C, size(descriptors,1), 1);
                    SST = vectDiff' * vectDiff;
                    SST = trace(SST);
      %             SST = SSW + SSB;

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
                   CInter = (sum(sum(pdist(centers)))*2) / (clusterCount^2);

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

%                    % Finally, calculate DB index.
%                    DB = 0;
%                    for centerItr = 1:clusterCount
%                         tempArr = zeros(clusterCount,1);
%                         for centerItr2 = 1:clusterCount
%                              if centerItr == centerItr2
%                                   continue;
%                              end
%                              tempArr(centerItr2) = (sigmas(centerItr) + sigmas(centerItr2)) / sqrt(sum((centers(centerItr,:) - centers(centerItr2,:)).^2));
%                         end
%                         DB = DB + max(tempArr);
%                    end
%                    DB = DB / clusterCount;
%                    DBVals(stepItr) = DB;
                   stepItr = stepItr + 1;
              end

              % Show DB index.
              figure, plot(sampleCounts, cutoffRatios);
              saveas(gcf, [options.debugFolder '/level' num2str(levelItr) '_ValidityIdx.png']);
              close all;
%               figure, plot(sampleCounts, DBVals);
%               saveas(gcf, [options.debugFolder '/level' num2str(levelItr) '_DBIdx_NoNorm.png']);
%               close all;
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
                   clusterCount = options.reconstruction.numberOfORNodes;
     %              [~, idx] = max(dunnVals);
     %              clusterCount = sampleCounts(idx(1));
              end
         else
              clusterCount = options.reconstruction.numberOfORNodes;
         end
         
         % Find optimal number of clusters.
%          [~, maxIdx] = max(dunnVals);
%          clusterCount = sampleCounts(maxIdx);

         display(['........ Obtained ' num2str(clusterCount) ' clusters! Finishing clustering.']);
         clusters = cluster(Z, 'maxclust', clusterCount);

         % Visualize dendogram.
          try %#ok<TRYNC>
               figure, dendrogram(Z, min(clusterCount, 50));
               saveas(gcf, [options.debugFolder '/level' num2str(levelItr) '_dendogram.png']);
          end
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
    
    clearvars -except vocabLevel graphLevel newDistanceMatrix
end