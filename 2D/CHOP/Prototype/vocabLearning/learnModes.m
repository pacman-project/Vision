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
%>               modes are of the form: [ node11, node12, coord11, coord12;
%>                                        node21, node22, coord21, coord22;
%>                                      ...]
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 21.01.2014
function [modes] = learnModes(currentLevel, newDistanceMatrix, edgeCoords, edgeIdMatrix, datasetName, levelItr, currentFolder)
    display('Learning modes...');
    maxSamplesPerMode = 500;
    minSamplesPerMode = 10;   
    maximumModes = 50;
    dummySigma = 0.1;
    halfSize = ceil(size(edgeIdMatrix,1) / 2);
    
    % Set initial data structures for processing 
    edges = cat(1, currentLevel.adjInfo);
    nodeIds = [currentLevel.labelId]';
    nodeCoords = cat(1, currentLevel.precisePosition);
    
    % Create redundancy vector, matching or nodes.
    if ~isempty(newDistanceMatrix)
         IC = 1:size(newDistanceMatrix,1);
         for clusterItr = 1:size(newDistanceMatrix,1)
              assignedCluster = find(newDistanceMatrix(clusterItr, 1:clusterItr) < 0.001, 1, 'first');
              if assignedCluster ~= clusterItr
                   IC(clusterItr) = assignedCluster;
              end
         end
    else
         IC = 1:numel(unique(nodeIds));
    end
    
    % If no edges exist, return.
    if isempty(edges)
       modes = [];
       return; 
    end
    
    % Anonymize edges, we're interested in labels, not ids.
    allEdges = [IC(nodeIds(edges(:,1:2))), int32(edgeCoords(edges(:,3), :))];
    
%    % Eliminate edges for which second node's labels are smaller than first.
%     allEdges = allEdges(allEdges(:,1) <= allEdges(:,2),:);
    
    %% Learn unique edge types and put them in cells for fast processing using parfor.
    [uniqueEdgeTypes, ~, IA] = unique(allEdges(:,1:2), 'rows');
    numberOfUniqueEdges = size(uniqueEdgeTypes,1);
    modes = cell(numberOfUniqueEdges,1);
    uniqueEdgeSamples = cell(numberOfUniqueEdges,1);
    uniqueEdgeCoords = cell(numberOfUniqueEdges,1);
    parfor uniqueEdgeItr = 1:numberOfUniqueEdges
        samplesForEdge = allEdges(IA==uniqueEdgeItr,3:4); %#ok<PFBNS>
        sampleIds = edges(IA==uniqueEdgeItr,1:2);
        sampleEdgeCoords = nodeCoords(sampleIds(:,2),:) - nodeCoords(sampleIds(:,1),:);
        
        %% If there are too many samples, get random samples.
        if size(samplesForEdge,1)>maxSamplesPerMode
            [samplesForEdge, idx] = datasample(samplesForEdge, maxSamplesPerMode, 'Replace', false);
            sampleEdgeCoords = sampleEdgeCoords(idx,:);
        end
        uniqueEdgeSamples(uniqueEdgeItr) = {samplesForEdge};
        uniqueEdgeCoords(uniqueEdgeItr) = {sampleEdgeCoords};
    end
    
    %% For each unique edge type (node1-node2 pair), estimate modes and save them in modes array.
    parfor uniqueEdgeItr = 1:numberOfUniqueEdges
  %      display(num2str(uniqueEdgeItr));
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
        statistics = zeros(numberOfClusters, 11, 'single');
        statistics(:,1:2) = single(repmat(edgeType, numberOfClusters, 1));
       for centerItr = 1:numberOfClusters
          clusterSamples = samples(classes==centerItr,:);
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
       
        %% We should also visualize the Gaussian PDFs.
        [X1, X2] = meshgrid(1:size(edgeIdMatrix,1), 1:size(edgeIdMatrix,1));
        X1 = (size(edgeIdMatrix,1) + 1) - X1;
        F = zeros(numel(X1), 1, class(statistics));
        gaussImg = zeros(size(edgeIdMatrix));
        for centerItr = 1:numberOfClusters
             try
                 curModePoints = mvnpdf([X1(:)-halfSize, X2(:)-halfSize],statistics(centerItr,4:5) , [statistics(centerItr,6:7); statistics(centerItr,8:9)]);
             catch
                 curModePoints = mvnpdf([X1(:)-halfSize, X2(:)-halfSize],statistics(centerItr,4:5) , statistics(centerItr,[6,9]));
             end
             curModePoints = curModePoints/max(max(curModePoints));
             F = max(F, curModePoints);
        end
        gaussImg(:) = F;
        gaussImg = imrotate(gaussImg, 90);
        gaussImg = uint8(round(gaussImg * 255));
        gaussImg(gaussImg==0) = 1;
        gaussImg = label2rgb(gaussImg, 'jet');
        
          %% Print the modes.
          distributionImg = zeros(size(edgeIdMatrix));
          samplesToWrite = floor(samples + halfSize) ;

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
%    modes = cat(1, modes{:});
    numberOfVocabNodes = numel(IC);
    newModes = cell(numberOfVocabNodes * numberOfVocabNodes, 1);
    for clusterItr = 1:numberOfVocabNodes
         for clusterItr2 = 1:numberOfVocabNodes
              idx = find(uniqueEdgeTypes(:,1) == IC(clusterItr) & uniqueEdgeTypes(:,2) == IC(clusterItr2));
              if ~isempty(idx)
                   idx = idx(1);
                   tempModes = modes{idx};
                   tempModes(:,1) = clusterItr;
                   tempModes(:,2) = clusterItr2;
                   newModes{(clusterItr-1) * numberOfVocabNodes + clusterItr2} = tempModes;
              end
         end
    end
    modes = cat(1, newModes{:});
        
%     %% Add reverse modes to the modes array.
%     reversedModes = modes(modes(:,1) ~= modes(:,2),:);
%     tempArr = reversedModes(:,1);
%     reversedModes(:,1) = reversedModes(:,2);
%     reversedModes(:,2) = tempArr;
%     reversedModes(:,3) = (double(max(max(edgeIdMatrix))) - reversedModes(:,3)) + 1;
%     reversedModes(:,4:5) = reversedModes(:,4:5) * -1;
%     modes = [modes; reversedModes];

    % Sort array.
    modes = sortrows(modes);
end