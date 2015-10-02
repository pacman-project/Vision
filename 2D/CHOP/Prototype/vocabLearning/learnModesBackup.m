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
function [modes] = learnModes(currentLevel, edgeCoords, edgeIdMatrix, datasetName, levelItr, currentFolder)
    display('Learning modes...');
    maxSamplesPerMode = 200;
    minSamplesPerMode = 2;   
    maximumModes = 4;
    dummySigma = 0.1;
    halfSize = ceil(size(edgeIdMatrix,1) / 2);
    edgeQuantize = size(edgeIdMatrix,1);
    
    % Set initial data structures for processing 
    edges = cat(1, currentLevel.adjInfo);
    nodeIds = [currentLevel.labelId]';
    nodeCoords = cat(1, currentLevel.precisePosition);
    
    % If no edges exist, return.
    if isempty(edges)
       modes = [];
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
    for uniqueEdgeItr = 1:numberOfUniqueEdges
        samplesForEdge = allEdges(IA==uniqueEdgeItr,3:4); %#ok<PFBNS>
        sampleIds = edges(IA==uniqueEdgeItr,1:2); %#ok<PFBNS>
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
    for uniqueEdgeItr = 1:numberOfUniqueEdges
        w = warning('off', 'all');
        samples = single(uniqueEdgeSamples{uniqueEdgeItr});
        sampleCoords = single(uniqueEdgeCoords{uniqueEdgeItr});
        edgeType = uniqueEdgeTypes(uniqueEdgeItr,:);
        
        %% Assign a label to each sample.
        if maximumModes == 1
            classes = ones(size(samples,1),1);
        else
            classes = assignModes(samples, minSamplesPerMode, maximumModes);
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
          normalizedCenter = round(statistics(centerItr,4:5)) + halfSize;
          statistics(centerItr,3) = double(edgeIdMatrix(normalizedCenter(1), normalizedCenter(2))); %#ok<PFBNS>
          
          if size(clusterSamples,1) > 1
               covMat = cov(clusterSamples);
               statistics(centerItr,6:7) = covMat(1,:);
               statistics(centerItr,8:9) = covMat(2,:);
               if nnz(statistics(centerItr,6:9)) < 2
                    statistics(centerItr,[6, 9]) = max(statistics(centerItr,[6, 9]), [dummySigma, dummySigma]);
               end
          else
               statistics(centerItr,[6, 9]) = dummySigma;
          end
       end
       clusterWeights = numberOfSamplesPerCluster / sum(numberOfSamplesPerCluster);
       statistics(:,12) = clusterWeights;
        
        %% Here, we create probabilities of each entry based on each mode of the distribution.
        % First, we obtain all possible points in the receptive field.
         yArr = repmat(1:edgeQuantize, edgeQuantize, 1);
         xArr = repmat((1:edgeQuantize)', 1, edgeQuantize);
        points = [xArr(:), yArr(:)];
        
        % In order to calculate probabilities, we define an area for each
        % point. 
        points = (points-(1/2)) / edgeQuantize;
        endPoints = points + 1/edgeQuantize;
        topPoints = points;
        topPoints(:,1) = points(:,1) + 1/edgeQuantize;
        rightPoints = points;
        rightPoints(:,2) = points(:,2) + 1/edgeQuantize;
        
        % Obtain probabilities based on each mode.
        allProbs = zeros(numberOfClusters, edgeQuantize, edgeQuantize, 'single');
        for centerItr = 1:numberOfClusters
            try
                probs = mvncdf(points, statistics(centerItr,4:5) , [statistics(centerItr,6:7); statistics(centerItr,8:9)]);
            catch %#ok<*CTCH>
                probs = mvncdf(points, statistics(centerItr,4:5) , statistics(centerItr,[6,9]));
            end
            try
                endProbs = mvncdf(endPoints, statistics(centerItr,4:5) , [statistics(centerItr,6:7); statistics(centerItr,8:9)]);
            catch
                endProbs = mvncdf(endPoints, statistics(centerItr,4:5) , statistics(centerItr,[6,9]));
            end
            try
                topProbs = mvncdf(topPoints, statistics(centerItr,4:5) , [statistics(centerItr,6:7); statistics(centerItr,8:9)]);
            catch
                topProbs = mvncdf(topPoints, statistics(centerItr,4:5) , statistics(centerItr,[6,9]));
            end 
            try
                rightProbs = mvncdf(rightPoints, statistics(centerItr,4:5) , [statistics(centerItr,6:7); statistics(centerItr,8:9)]);
            catch
                rightProbs = mvncdf(rightPoints, statistics(centerItr,4:5) , statistics(centerItr,[6,9]));
            end
            
            % Calculate point probs based on the integral.
            pointProbs = (endProbs - (topProbs + rightProbs)) + probs;
            pointProbs = pointProbs * statistics(centerItr,12);
            allProbs(centerItr, :,:) = reshape(pointProbs, size(edgeIdMatrix));
        end
        
        % Get the probability image for visualization.
        gaussImg = squeeze(max(allProbs, [], 1));
        gaussImg = gaussImg / max(max(gaussImg));
        gaussImg = imrotate(gaussImg, 90);
        gaussImg = uint8(round(gaussImg * 255));
        gaussImg(gaussImg==0) = 1;
        gaussImg = label2rgb(gaussImg, 'jet');
        
          %% Print the modes.
          distributionImg = zeros(size(edgeIdMatrix));
          samplesToWrite = floor(samples * edgeQuantize) ;

          % If no samples are to be written, move on.E
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
          
        %% Save stats and move on.
        modes(uniqueEdgeItr) = {statistics};
        warning(w);
    end
    clear uniqueEdgeSamples;
    modes = cat(1, modes{:});

    % Sort array.
    modes = sortrows(modes);
end