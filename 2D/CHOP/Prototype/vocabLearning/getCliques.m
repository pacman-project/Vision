function [jointPositions, savedCliques, uniqueSubParts, numberOfSubParts] = getCliques(candidateParts, allCliques, firstInstances, prevLevelPositions, allInstancePositions, rfSize, partItr)
  validSubPartCount = nnz(candidateParts(partItr,:));
  relevantCliques = allCliques(firstInstances(partItr):(firstInstances(partItr+1)-1),:);
  instancePositions = allInstancePositions(firstInstances(partItr):(firstInstances(partItr+1)-1),:);
  relevantCliques = relevantCliques(:, (end + (1 - validSubPartCount)):end);
  maxPointsToCluster = 200;

  %% Learn feature maps here. Statistical maps are generated from instances.
  % The feature map has as many layers as the number of individual
  % sub labels in children.
  subParts = candidateParts(partItr,:);
  subParts = subParts(subParts>0);
  numberOfSubParts = numel(subParts);
  uniqueSubParts = fastsortedunique(subParts);

  % In order to learn joint probabilities, we need to re-order the
  % nodes if there are duplicate children ( Children having the
  % same label). In order to do this, we cluster children in
  % duplicate channels, and re-order them after assigning every
  % children to a class. Children of the same class are ordered in
  % the same way. Outliers are removed.
  if numel(uniqueSubParts) ~= numberOfSubParts
       % Learn each feature map by itself.
       for mapItr = 1:numel(uniqueSubParts)
            % If this type of object has been seen only once, we
            % don't need to process anything.
            if nnz(subParts == uniqueSubParts(mapItr)) == 1
                 continue;
            end

            % Get relevant instances.
            relevantRows = subParts == uniqueSubParts(mapItr);
            allRelevantChildren = relevantCliques(:,relevantRows);

            % Select a subset of the data to cluster.
            if size(allRelevantChildren,1) > maxPointsToCluster 
                 selectedRows = datasample(1:size(allRelevantChildren,1), maxPointsToCluster, 'Replace', false);
            else
                 selectedRows = 1:size(allRelevantChildren,1);
            end

            % Collect data and perform clustering.
            allPositions = cell(nnz(relevantRows),1);
            for childItr = 1:nnz(relevantRows)
                 column = allRelevantChildren(:, childItr);
                 relativePositions = prevLevelPositions(column,:) - instancePositions;
                 allPositions{childItr} = relativePositions;
            end
            allPositions = cat(1, allPositions{:});

            % Perform minimum entropy-based clustering.
            clusteredPositions = allPositions(selectedRows,:);
%                     [~, clusterCenters] = gmeans(double(clusteredPositions), maxClusters, 0.01);
%                     
%                     % Plot the points in space and write stuff down.
%                     [~, ~, repIdx] = unique(clusteredPositions, 'rows');
%                     repetitions = hist(repIdx, unique(repIdx));
%                     repetitions = repetitions/max(repetitions);
%                     figure, plot3(clusteredPositions(:,1), clusteredPositions(:,2), repetitions(repIdx), 'rx');

%                    [~, clusterCenters2] = gmeans(double(allPositions), maxClusters, 0.0001);
            if size(clusteredPositions,1) == 1
                 clusters = 1;
            else
                 [clusters, ~] = gmeans(double(clusteredPositions*(1/((rfSize-1) * sqrt(size(clusteredPositions,2))))), nnz(subParts == uniqueSubParts(mapItr)), 0.05, @checkGaussian);
            end
            
%                    clusters = mec(clusteredPositions*(2/realRFSize), 'c', numberOfSubParts, 'kmeans_i', 5, 'kernel', 1, 'mec_maxiter', 20, 'width', 3);
            numberOfClusters  = numel(unique(clusters));
            clusterCenters = zeros(numberOfClusters,2, 'single');
            for clusterItr = 1:numberOfClusters
                 clusterCenters(clusterItr,:) = sum(clusteredPositions(clusters == clusterItr,:),1) / nnz(clusters == clusterItr);
            end

            % Assign all existing samples to clusters
            distances = pdist2(allPositions, single(clusterCenters));
            [~,allClusters] = min(distances, [], 2);
            allClusters = reshape(allClusters, [size(relevantCliques,1), nnz(relevantRows)]);
            [~, sortIdx] = sort(allClusters,2);

            % Sort the rows and write them back.
            for instanceItr = 1:size(allRelevantChildren,1)
                 allRelevantChildren(instanceItr,:) = allRelevantChildren(instanceItr, sortIdx(instanceItr,:));
            end
            relevantCliques(:,relevantRows) = allRelevantChildren;
       end
  end
  % Calculate joint positions.
  jointPositions = cell(numberOfSubParts,1);
  for childItr = 1:numberOfSubParts
       column = relevantCliques(:, childItr);
       relativePositions = prevLevelPositions(column,:) - instancePositions;
       jointPositions{childItr} = relativePositions;
  end
  jointPositions = cat(2, jointPositions{:});
  
  % Save cliques along with their padding.
  savedCliques = cat(2, zeros(size(relevantCliques,1), size(allCliques,2)-size(relevantCliques,2)), relevantCliques);
end

