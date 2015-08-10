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
function [modes] = learnModes(currentLevel, edgeCoords, edgeIdMatrix)
    display('Learning modes...');
    maxSamplesPerMode = 200;
    minSamplesPerMode = 4;   
    maximumModes = 50;
    dummySigma = 0.1;
    halfSize = ceil(size(edgeIdMatrix,1) / 2);
    
    % Set initial data structures for processing 
    edges = cat(1, currentLevel.adjInfo);
    nodeIds = [currentLevel.labelId]';
    
    % If no edges exist, return.
    if isempty(edges)
       modes = [];
       return; 
    end
    
    % Anonymize edges, we're interested in labels, not ids.
    allEdges = [nodeIds(edges(:,1:2)), int32(edgeCoords(edges(:,3), :))];
    
%    % Eliminate edges for which second node's labels are smaller than first.
%     allEdges = allEdges(allEdges(:,1) <= allEdges(:,2),:);
    
    %% Learn unique edge types and put them in cells for fast processing using parfor.
    [uniqueEdgeTypes, ~, IA] = unique(allEdges(:,1:2), 'rows');
    numberOfUniqueEdges = size(uniqueEdgeTypes,1);
    modes = cell(numberOfUniqueEdges,1);
    uniqueEdgeSamples = cell(numberOfUniqueEdges,1);
    parfor uniqueEdgeItr = 1:numberOfUniqueEdges
        samplesForEdge = allEdges(IA==uniqueEdgeItr,3:4); %#ok<PFBNS>
        %% If there are too many samples, get random samples.
        if size(samplesForEdge,1)>maxSamplesPerMode
            samplesForEdge = datasample(samplesForEdge, maxSamplesPerMode, 'Replace', false);
        end
        uniqueEdgeSamples(uniqueEdgeItr) = {samplesForEdge};
    end
    
    %% For each unique edge type (node1-node2 pair), estimate modes and save them in modes array.
    parfor uniqueEdgeItr = 1:numberOfUniqueEdges
  %      display(num2str(uniqueEdgeItr));
        w = warning('off', 'all');
        samples = double(uniqueEdgeSamples{uniqueEdgeItr});
        edgeType = uniqueEdgeTypes(uniqueEdgeItr,:);
        
        %% Assign a label to each sample.
        if maximumModes == 1
            classes = ones(size(samples,1),1);
        else
            classes = assignModes(samples, minSamplesPerMode, maximumModes);
        end

        %% Calculate statistics (mu, sigma)
        numberOfClusters = max(classes);
        statistics = zeros(numberOfClusters, 9);
        statistics(:,1:2) = double(repmat(edgeType, numberOfClusters, 1));
       for centerItr = 1:numberOfClusters
          clusterSamples = samples(classes==centerItr,:);
          statistics(centerItr,4:5) = mean(clusterSamples,1);
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
        
        modes(uniqueEdgeItr) = {statistics};
           
        warning(w);
    end
    clear uniqueEdgeSamples;
    modes = cat(1, modes{:});
        
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