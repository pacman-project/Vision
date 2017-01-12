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
function [vocabLevel, vocabularyDistributions, graphLevel, newDistanceMatrix] = postProcessParts(vocabLevel, graphLevel, vocabularyDistributions, levelItr, options)
    distType = options.distType;
    stopFlag = options.stopFlag;
    smallImageSize = 50;
    % Get filter size.
    filterSize = size(options.filters{1});
    
    %% As a case study, we replace gabors with 1D gaussian filters    
    visFilters = options.filters;
    visFilters = cellfun(@(x) (x - min(min(x))) / (max(max(x)) - min(min(x))), visFilters, 'UniformOutput', false);
    for filterItr = 1:numel(visFilters)
         visFilters{filterItr} = uint8(round(visFilters{filterItr} * 255));
    end
    visFilters = cat(3, visFilters{:});
    visFilters(visFilters<1) = 1;
    visFilters = double(visFilters);
    
    % Assign new labels of the remaining realizations.
    [remainingComps, ~, IC] = unique([graphLevel.labelId]);
    IC = num2cell(int32(IC));
    [graphLevel.labelId] = deal(IC{:});
    [graphLevel.realLabelId] = deal(IC{:});
    clear IC;
    
    % Eliminate unused compositions from vocabulary.
    vocabLevel = vocabLevel(remainingComps);
    vocabLevelDistributions = vocabularyDistributions{levelItr};
    vocabLevelDistributions = vocabLevelDistributions(remainingComps);
    vocabularyDistributions{levelItr} = vocabLevelDistributions;
    newLabelArr = num2cell(int32(1:numel(vocabLevel)));
    [vocabLevel.label] = deal(newLabelArr{:});
    
    %% Find the distance matrix among the remaining parts in vocabLevel.
   numberOfNodes = numel(vocabLevel);
   vocabChildren = {vocabLevel.children};
   largeSubIdx = cellfun(@(x) numel(x) > 1, vocabChildren);
   
   %% We  are experimenting with different distance functions.
   % Find the correct image size.
   imageSize = getRFSize(options, levelItr);
   imageSize = imageSize + filterSize - 1;
   imageCenter = floor(imageSize(1)/2)+1;
   
   %% First, for efficiency, we obtain pixel-level predictions for every part.
   level1Experts = cell(numberOfNodes, 1);
   for vocabNodeItr = 1:numberOfNodes
        % Backproject nodes using modal reconstructions.
        nodes = [vocabNodeItr, 0, 0, levelItr];
        experts = projectNode(nodes, vocabularyDistributions, 'modal', options);

        % Center the nodes.
        experts = double(experts);
        experts(:,2:3) = experts(:,2:3) + imageCenter;
        level1Experts{vocabNodeItr} = experts;
   end
   
   % Comparison of modal reconstructions involves creating a pixel
   % prediction for every pixel, and then looking for matches.
   muImgs = zeros(numberOfNodes, imageSize(1), imageSize(2));
   smallImageSize = min(smallImageSize, imageSize(1));
   smallMuImgs = zeros(numberOfNodes, smallImageSize, smallImageSize);
   newDistanceMatrix = zeros(numel(vocabLevel), 'single');

   % Create a blurring filter.
   H = fspecial('gaussian', 3, 1);
   
   % Get product of expert predictions.
   display('........ Visualizing parts using Product of Experts..');
   mkdir([options.debugFolder '/level' num2str(levelItr) '/modalProjection/']);
   debugFolder = options.debugFolder;
   parfor vocabNodeItr = 1:numberOfNodes
        [muImg, ~, ~] = obtainPoE(level1Experts{vocabNodeItr}, [], [], imageSize, visFilters, []);
        muImg = uint8(round(255*(double(muImg) / double(max(max(muImg))))));
        blurredMuImg = uint8(imfilter(muImg, H, 'replicate'));
        muImgs(vocabNodeItr,:,:) = blurredMuImg;
        
        % Save the image in a small array.
        if imageSize(1) > smallImageSize
            smallMuImgs(vocabNodeItr, :, :) = imresize(muImg, [smallImageSize, smallImageSize]); %#ok<PFOUS>
        else
            smallMuImgs(vocabNodeItr, :, :) = muImg;
        end
           
        blurredMuImg = uint8(round(255*(double(blurredMuImg) / double(max(max(blurredMuImg))))));
        imwrite(muImg, [debugFolder '/level' num2str(levelItr) '/modalProjection/' num2str(vocabNodeItr) '.png']);
        imwrite(blurredMuImg, [debugFolder '/level' num2str(levelItr) '/modalProjection/' num2str(vocabNodeItr) '_blurred.png']);
   end
   save([debugFolder '/level' num2str(levelItr) '/modalProjection/muImgs.mat'], 'smallMuImgs'); 
   
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
   elseif strcmp(distType, 'sift')
        descriptors = zeros(numel(vocabLevel), 128);
        tempSize = size(muImgs);
        firstSide = floor(tempSize(end-1) / 4);
        secSide = floor(tempSize(end)/4);
        imScale = min(firstSide, secSide);
        centerPoint = round([tempSize(end-1); tempSize(end)]/2);
        for vocabNodeItr = 1:numel(vocabLevel)
               [Ix, Iy] = vl_grad(im2double(squeeze(muImgs(vocabNodeItr,:,:))));
               mod      = sqrt(Ix.^2 + Iy.^2) ;
               ang      = atan2(Iy,Ix) ;
               grd      = shiftdim(cat(3,mod,ang),2) ;
               grd      = single(grd) ;
               descriptors(vocabNodeItr,:) = single(vl_siftdescriptor(grd, [centerPoint; 1; 0], 'Magnif', imScale)')/255;
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
   elseif strcmp(distType, 'context')
        % Shape context matching.
        % We extract a number of points (50) from every object contour, and then
        % match them. 
        defaultoptions=struct('r_max',4,'r_min',1e-3,'r_bins',15,'a_bins',15,'rotate',0,'method',1,'maxdist',5);
        descriptors = zeros(numel(vocabLevel), defaultoptions.r_bins * defaultoptions.a_bins);
        for vocabNodeItr = 1:numel(vocabLevel)
             img = squeeze(muImgs(vocabNodeItr,:,:));
             binaryImg = img > max(max(img))/2;
             points = find(binaryImg);
             if numel(points) > 100
                  points = points(randperm(size(points,1), 100));
             end
             [pointsX, pointsY] = ind2sub(size(img), points);
             allPoints = [pointsX, pointsY; 0.5+size(img)/2];
             
             % Create a shape context descriptor for this part.
             [Points1,~]=NormalizePoints(allPoints,allPoints,defaultoptions);
             F1=getHistogramFeatures(Points1(end,:),Points1(1:(end-1),:),[],defaultoptions);
             descriptors(vocabNodeItr,:) = F1';
        end
        newDistanceMatrix = pdist3(descriptors, descriptors, 'chisq');
   end
   
   % Renormalize distance matrix.
   if nnz(newDistanceMatrix) > 0
       newDistanceMatrix = single(newDistanceMatrix/max(max(newDistanceMatrix)));
   end
   if isempty(newDistanceMatrix)
        newDistanceMatrix = 0;
   end
    
    %% Find an optimal number of clusters based on Davies-Bouldin index.
    if size(newDistanceMatrix,1) < 3
         clusters = (1:size(newDistanceMatrix,1))';
    else
         % Convert the distance matrix to vector form.
         newDistanceMatrixValid = newDistanceMatrix(largeSubIdx, largeSubIdx);
         newDistanceMatrixVect = squareform(newDistanceMatrixValid);
         
         %% Finally, we implement the OR nodes here.
         % We apply agglomerative clustering on the generated distance matrix,
         % and group parts based on their similarities. We have a limited number
         % of resources when selecting parts.
         % First, we check for the necessity.
         % All good, group the nodes here.
         if nnz(largeSubIdx) > 1
            Z = linkage(newDistanceMatrixVect, 'ward');
         end
         
         % Assign a fix value for the class count.
         if ~stopFlag
             clusterCount = min(numberOfNodes, max(options.reconstruction.minNumberOfORNodes, round(nnz(largeSubIdx)/options.partCountDenom)));
%             clusterCount = numberOfNodes;
         else
             clusterCount = numberOfNodes;
         end

         display(['........ Obtained ' num2str(clusterCount) ' clusters! Finishing clustering.']);
         if nnz(largeSubIdx) > 1
            clusters = cluster(Z, 'maxclust', clusterCount);
         else
            clusters = 1; 
         end
         
         % Finally, add single node subs.
         finalClusters = zeros(numberOfNodes, 1);
         finalClusters(largeSubIdx) = clusters;
         finalClusters(~largeSubIdx) = (max(clusters) + 1) + (1:(numel(largeSubIdx) - nnz(largeSubIdx)))';
         clusters = finalClusters; 

         % Visualize dendogram.
          try %#ok<TRYNC>
              if usejava('jvm')
                  figure, dendrogram(Z, min(clusterCount, 50));
                  saveas(gcf, [options.debugFolder '/level' num2str(levelItr) '_dendogram.png']);
              end
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
    
    clearvars -except vocabLevel vocabularyDistributions graphLevel newDistanceMatrix
end