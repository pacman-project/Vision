%> Name: learnModes
%>
%> Description: Given the node list of all images, this function learns the
%> modes by clustering pairwise relative positions in 2D space. 
%>
%> @param mainGraph The object graphs' data structure.
%> @param options Program options.
%> @param currentLevel The current scene graph level.
%> @param datasetName Name of the dataset.
%> 
%> @retval modes The mode list representing edge categories.
%>               modes are of the form: [ nodeLabel1, nodeLabel2, edgeId, rfCoord1, rfCoord2, cov11, cov12, cov21, cov22, coord1, coord2, weight;
%>                                      ...]
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 21.01.2014
function [modes, modeProbArr] = learnModes(currentLevel, edgeCoords, edgeIdMatrix, datasetName, levelItr, currentFolder, debug)
    display('Learning modes...');
    maxSamplesPerMode = 200;
    minSamplesPerMode = 2;   
    if levelItr == 1
         maximumModes = 10;
    else
         maximumModes = 5;
    end
    dummySigma = 0.02;
    noiseSigma = 0.001;
    minProb = 0.001;
    halfSize = ceil(size(edgeIdMatrix,1) / 2);
    edgeQuantize = size(edgeIdMatrix,1);
    
    cx=halfSize;cy=cx;ix=edgeQuantize;iy=ix;r=halfSize;
    [x,y]=meshgrid(-(cx-1):(ix-cx),-(cy-1):(iy-cy));
    c_mask=((x.^2+y.^2)<r^2);
    
    % Set initial data structures for processing 
    edges = cat(1, currentLevel.adjInfo);
    nodeIds = [currentLevel.labelId];
    nodeCoords = cat(1, currentLevel.precisePosition);
    
    % If no edges exist, return.
    if isempty(edges)
       modes = [];
       modeProbArr = [];
       return; 
    end
    
    % Anonymize edges, we're interested in labels, not ids.
    allEdges = [nodeIds(edges(:,1:2)), int32(edgeCoords(edges(:,3), :))];
    
    %% Learn unique edge types and put them in cells for fast processing using parfor.
    [uniqueEdgeTypes, ~, IA] = unique(allEdges(:,1:2), 'rows');
    numberOfUniqueEdges = size(uniqueEdgeTypes,1);
    modes = cell(numberOfUniqueEdges,1);
    uniqueEdgeSamples = cell(numberOfUniqueEdges,1);
    uniqueEdgeCoords = cell(numberOfUniqueEdges,1);
    parfor uniqueEdgeItr = 1:numberOfUniqueEdges
        tempIdx = IA==uniqueEdgeItr;
        samplesForEdge = allEdges(tempIdx,3:4); %#ok<PFBNS>
        sampleIds = edges(tempIdx,1:2); %#ok<PFBNS>
        sampleEdgeCoords = nodeCoords(sampleIds(:,2),:) - nodeCoords(sampleIds(:,1),:); %#ok<PFBNS>
        
        %% If there are too many samples, get random samples.
        if size(samplesForEdge,1)>maxSamplesPerMode
            [samplesForEdge, idx] = datasample(samplesForEdge, maxSamplesPerMode, 'Replace', false);
            sampleEdgeCoords = sampleEdgeCoords(idx,:);
        end
        %Downsample samplesForEdge and save them.
        samplesForEdge = (double(samplesForEdge + halfSize) / edgeQuantize);
        uniqueEdgeSamples(uniqueEdgeItr) = {samplesForEdge};
        uniqueEdgeCoords(uniqueEdgeItr) = {sampleEdgeCoords};
    end
    
    %% For each unique edge type (node1-node2 pair), estimate modes and save them in modes array.
    modeProbArr = cell(numberOfUniqueEdges, 1);
    for uniqueEdgeItr = 1:numberOfUniqueEdges
        probs = [];endProbs=[];topProbs=[];rightProbs=[]; %#ok<NASGU>
        w = warning('off', 'all');
        samples = single(uniqueEdgeSamples{uniqueEdgeItr});
        sampleCoords = single(uniqueEdgeCoords{uniqueEdgeItr});
        edgeType = uniqueEdgeTypes(uniqueEdgeItr,:);
        
        %% Assign a label to each sample.
        if maximumModes == 1
            classes = ones(size(samples,1),1);
        else
            classes = assignModes(samples, minSamplesPerMode, maximumModes);
            tempNumberOfClusters = max(classes);
            if tempNumberOfClusters > 1
                 classes = kmeans(samples, tempNumberOfClusters, 'EmptyAction', 'singleton');
            end
        end

        %% Calculate statistics (mu, sigma)
        numberOfClusters = max(classes);
        numberOfSamplesPerCluster = zeros(numberOfClusters,1);
        statistics = zeros(numberOfClusters, 12, 'single');
        statistics(:,1:2) = single(repmat(edgeType, numberOfClusters, 1));
       for centerItr = 1:numberOfClusters
          clusterSamples = samples(classes==centerItr,:);
          numberOfSamplesPerCluster(centerItr) = size(clusterSamples,1);
          clusterCoords = sampleCoords(classes == centerItr, :);
          statistics(centerItr,4:5) = mean(clusterSamples,1);
          statistics(centerItr,10:11) = mean(clusterCoords,1);
          normalizedCenter = round(statistics(centerItr,4:5) * edgeQuantize);
          statistics(centerItr,3) = double(edgeIdMatrix(normalizedCenter(1), normalizedCenter(2))); %#ok<PFBNS>
          
          if size(clusterSamples,1) > 1
               %% Calculate covariance matrix.
               covMat = cov(clusterSamples);
               covMat([1,4]) = max(covMat([1,4]), dummySigma);
 %              minChannelVar = max(minVarMultiplier)
               
               % Smooth the covariance matrix, in case we're too narrow on
               % one dimension.
               covMat = max(covMat, 0);
               try
                    dummyProb = mvnpdf(statistics(centerItr,4:5), statistics(centerItr,4:5), covMat); %#ok<NASGU>
               catch %#ok<CTCH>
                    covMat = covMat .* eye(size(covMat,1)) * noiseSigma;
               end
               
               % Assign matrix.
               statistics(centerItr,6:7) = covMat(1,:);
               statistics(centerItr,8:9) = covMat(2,:);
          else
               statistics(centerItr,[6, 9]) = dummySigma;
          end
       end
       clusterWeights = numberOfSamplesPerCluster / sum(numberOfSamplesPerCluster);
       statistics(:,12) = clusterWeights;
       
       % If there are two or more clusters on the same mean, keep only one mode,
       % discard the other(s).
       if numel(unique(statistics(:,3))) ~= numel(statistics(:,3))
            % Find valid clusters and cancel the rest.
            [~, validClusters, ~] = unique(statistics(:,3), 'stable');
            numberOfClusters = numel(validClusters);
            statistics = statistics(validClusters,:);
            validSampleIdx = ismember(classes, validClusters);
            
            % Update samples and their classes.
            samples = samples(validSampleIdx,:);
            classes = classes(validSampleIdx);
            labelAssgnArr = zeros(max(validClusters),1, 'int32');
            labelAssgnArr(validClusters) = 1:numel(validClusters);
            classes = labelAssgnArr(classes);
            
            % Renormalize cluster weights.
            statistics(:,12) = statistics(:,12)/sum(statistics(:,12));
       end
        
        %% Here, we create probabilities of each entry based on each mode of the distribution.
        % First, we obtain all possible points in the receptive field.
        yArr = repmat(1:edgeQuantize, edgeQuantize, 1);
        xArr = repmat((1:edgeQuantize)', 1, edgeQuantize);
        points = [xArr(:), yArr(:)];
        
        % In order to calculate probabilities, we define an area for each
        % point. 
        points = points / edgeQuantize;
        numberOfPoints = size(points,1);
        halfPixelSize = 1/(2*edgeQuantize);
        halfPixelArr = ones(numberOfPoints,1) * halfPixelSize;
        
        % Obtain probabilities based on each mode.
        allProbs = zeros(numberOfClusters, edgeQuantize, edgeQuantize, 'single');
        for centerItr = 1:numberOfClusters
            % We use PDFs to calculate pseudo-statistics (faster, accurate
            % enough).
%            try
            lowProbs = mvncdf(points + [-halfPixelArr, -halfPixelArr], statistics(centerItr,4:5) , [statistics(centerItr,6:7); statistics(centerItr,8:9)]);
            leftProbs = mvncdf(points + [halfPixelArr, -halfPixelArr], statistics(centerItr,4:5) , [statistics(centerItr,6:7); statistics(centerItr,8:9)]);
            rightProbs = mvncdf(points + [-halfPixelArr, halfPixelArr], statistics(centerItr,4:5) , [statistics(centerItr,6:7); statistics(centerItr,8:9)]);
            highProbs = mvncdf(points + [halfPixelArr, halfPixelArr], statistics(centerItr,4:5) , [statistics(centerItr,6:7); statistics(centerItr,8:9)]);
            pointProbs = (highProbs - (leftProbs + rightProbs)) + lowProbs;
%            pointProbs = mvnpdf(points, statistics(centerItr,4:5) , [statistics(centerItr,6:7); statistics(centerItr,8:9)]);
%            catch %#ok<CTCH>
%                pointProbs = mvnpdf(points, statistics(centerItr,4:5) , statistics(centerItr,[6,9]));
%            end
            
            % Weight probabilities with the number of samples in each
            % cluster.
            pointProbs = pointProbs * statistics(centerItr,12);
            
            % Write probabilities into a matrix.
            allProbs(centerItr, :,:) = reshape(pointProbs, size(edgeIdMatrix));
        end
        initialAllProbs = allProbs;
        [probVals, idx] = max(allProbs, [], 1);
        clusterAreas = squeeze(idx);
        clusterAreas(~c_mask) = 0;
        clusterAreas(probVals < minProb) = 0;     
        
        % Assign minimum probabilities to the areas not covered by each
        % cluster.
        for centerItr = 1:numberOfClusters
             probMatrix = squeeze(allProbs(centerItr,:,:));
             probMatrix(clusterAreas ~= centerItr) = 0;
             
             % Obtain 
%             probMatrix(c_mask) = max(minProb, probMatrix(c_mask));
             allProbs(centerItr,:,:) = probMatrix;
        end
        
%        fillFlag = false;
%        contFlag = true;
%        iterations = 0;
%         while contFlag && iterations < numberOfClusters
%              contFlag = false;
%              % Update cluster areas so that all areas consist of a single
%              % connected component.
%              for centerItr = 1:numberOfClusters
%                   CC = bwconncomp(clusterAreas == centerItr, 4);
% 
%                   % If there's more than one 
%                   if CC.NumObjects> 1
%                        pixelCounts = cellfun(@(x) numel(x), CC.PixelIdxList);
%                        [~,selectedArea] = max(pixelCounts);
%                        otherPixList = cat(1, CC.PixelIdxList{setdiff(1:numel(pixelCounts), selectedArea)});
%                        fillFlag = true;
% 
%                        % Also mark probabilities as zero.
%                        probMatrix = squeeze(allProbs(centerItr, :, :));
%                        probMatrix(otherPixList) = 0;
%                        allProbs(centerItr, :, :) = probMatrix;
%                   end
%              end
%              if fillFlag
%                   [combinedProbs, idx] = max(allProbs, [], 1);
%                   clusterAreas = squeeze(idx);
%                   clusterAreas(squeeze(combinedProbs) == 0) = 0;
%                   contFlag = true;
%              end
%              iterations = iterations + 1;
%         end

        for centerItr = 1:numberOfClusters
             CC = bwconncomp(clusterAreas == centerItr, 8);

             % If there's more than one 
             if CC.NumObjects> 1
                  pixelCounts = cellfun(@(x) numel(x), CC.PixelIdxList);
                  [~,selectedArea] = max(pixelCounts);
                  otherPixList = cat(1, CC.PixelIdxList{setdiff(1:numel(pixelCounts), selectedArea)});
                  clusterAreas(otherPixList) = 0;
             end
        end
        
        % Renormalize probabilities so they sum up to 1.
        for centerItr = 1:numberOfClusters
             probMatrix = squeeze(initialAllProbs(centerItr,:,:));
             probMatrix(clusterAreas ~= centerItr) = 0;
             
             % Obtain 
 %            probMatrix(c_mask) = max(minProb, probMatrix(c_mask));
             allProbs(centerItr,:,:) = probMatrix;
        end
        [combinedProbs, ~] = max(allProbs, [], 1);
        combinedProbs = combinedProbs / sum(sum(combinedProbs));
        
        % Assign probabilities.
        modeProbArr{uniqueEdgeItr} = uint8(clusterAreas);
        
        % Get the probability image for visualization.
        if debug
             gaussImg = squeeze(combinedProbs);
             gaussImg = gaussImg / max(max(gaussImg));
             gaussImg = uint8(round(gaussImg * 255));
             gaussImg(gaussImg==0) = 1;
             gaussImg = label2rgb(gaussImg, 'jet');
        
               %% Print the modes.
               distributionImg = zeros(size(edgeIdMatrix));
               samplesToWrite = floor(samples * edgeQuantize) ;

               % If no samples are to be written, move on.
               if numel(samplesToWrite) < 1
                    continue;
               end
              samplesInd = sub2ind(size(distributionImg), samplesToWrite(:,1), samplesToWrite(:,2));
              distributionImg(samplesInd) = classes;

               if ~exist([currentFolder '/debug/' datasetName '/level' num2str(levelItr) '/pairwise/'], 'dir')
                    mkdir([currentFolder '/debug/' datasetName '/level' num2str(levelItr) '/pairwise/']);
               end
               sampleImg = label2rgb(distributionImg, 'jet', 'k', 'shuffle');
               imwrite(sampleImg, ...
                    [currentFolder '/debug/' datasetName '/level' num2str(levelItr) ...
                    '/pairwise/' num2str(edgeType(1)) '_' num2str(edgeType(2)) '_Samples.png']);
               imwrite(gaussImg, ...
                    [currentFolder '/debug/' datasetName '/level' num2str(levelItr) ...
                    '/pairwise/' num2str(edgeType(1)) '_' num2str(edgeType(2)) '_GaussMap.png']);
               areaImg = label2rgb(clusterAreas, 'jet', 'k', 'shuffle');
               imwrite(areaImg, ...
                    [currentFolder '/debug/' datasetName '/level' num2str(levelItr) ...
                    '/pairwise/' num2str(edgeType(1)) '_' num2str(edgeType(2)) '_ClusterAreas.png']);
        end
        %% Save stats and move on.
        modes(uniqueEdgeItr) = {statistics};
        warning(w);
    end
    modeProbArr = cat(3, modeProbArr{:});
    clear uniqueEdgeSamples;
    modes = cat(1, modes{:});

    % Sort array.
    modes = sortrows(modes);
    clearvars -except modes modeProbArr
end