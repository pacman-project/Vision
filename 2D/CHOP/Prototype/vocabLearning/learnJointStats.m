function [refinedClusterSamples, refinedCliques, refinedClusters, refinedInstancePositions, shownSamples, shownClusters, clusterDistributions, clusterThresholds, numberOfClusters, partClusters] = learnJointStats(jointPositions, instancePositions, savedCliques, rfSize, maxClusters) 
  % We will experiment with different algorithms here, 
  % from completely unsupervised to supervised.
  dummyStd = 0.0001;
  stdThr = 2.576;
  maxPointsToCluster = 500;
  if size(jointPositions,1) == 1
       partClusters = 1;
  else
       % Select a number of rows and cluster them, not all data.
       selectedRows = datasample(1:size(jointPositions,1), min(size(jointPositions,1), maxPointsToCluster), 'Replace', false);
       selectedPositions = jointPositions(selectedRows, :);
       [selectedClusters, ~] = gmeans(double(selectedPositions*(1/((rfSize-1) * sqrt(size(selectedPositions,2))))), maxClusters, 0.05, @checkGaussian);
       selectedCenters = zeros(max(selectedClusters), size(selectedPositions,2), 'single');
       for clusterItr = 1:max(selectedClusters)
            idx = selectedClusters == clusterItr;
            selectedCenters(clusterItr, :) = sum(selectedPositions(idx, :)) / nnz(idx);
       end
       distances = pdist2(jointPositions, selectedCenters);
       [~, partClusters] = min(distances, [], 2);
       if numel(unique(partClusters)) ~= max(partClusters)
            [~, ~, partClusters] = unique(partClusters);
       end
  end
  numberOfClusters  = numel(unique(partClusters));
  clusterCenters = zeros(numberOfClusters, size(jointPositions,2), 'single');
  clusterDistributions = cell(numberOfClusters,1);
  clusterThresholds = cell(numberOfClusters,1);
  refinedClusterSamples = cell(numberOfClusters,1);
  refinedInstancePositions = cell(numberOfClusters,1);
  refinedClusters = cell(numberOfClusters,1);
  refinedCliques = cell(numberOfClusters,1);
  eliminatedSamples = cell(numberOfClusters,1);
  validIdx = ones(numberOfClusters,1) > 0;
  savedClusterCount = 1;
  for clusterItr = 1:numberOfClusters
      idx = partClusters == clusterItr;
      clusterSamples = jointPositions(idx,:);
      clusterPositions = instancePositions(idx, :);
      clusterCliques = savedCliques(idx,:);
      clusterCenters(clusterItr,:) = sum(clusterSamples,1) / nnz(idx);

      % We eliminate nodes that are far away from the cluster centers
      % in every cluster.
      covMat = calculateCov(double(clusterSamples), dummyStd);
      distances = pdist2(clusterCenters(clusterItr,:), clusterSamples, 'mahalanobis', covMat);
      distances(isnan(distances)) = 0;
      thresh = max(min(3*stdThr+0.0001, min(distances)+0.0001), stdThr+0.0001); % 98 pct confidence interval
      validInstances = distances < thresh;
      
      if nnz(validInstances) == 0
           validIdx(clusterItr) = 0;
           continue;
      end

      % Filter out unlikely instances and save the data.
      refinedClusterSamples{clusterItr} = clusterSamples(validInstances,:);
      refinedCliques{clusterItr} = clusterCliques(validInstances, :);
      refinedInstancePositions{clusterItr} = clusterPositions(validInstances,:);
      eliminatedSamples{clusterItr} = clusterSamples(~validInstances, :);
      refinedClusters{clusterItr} = ones(nnz(validInstances),1) * savedClusterCount;

      % Save covariance matrix.
      clusterDistributions{clusterItr} = gmdistribution(clusterCenters(clusterItr,:), covMat);
      clusterThresholds{clusterItr} = double(thresh);
      
      % Save cluster count.
      savedClusterCount = savedClusterCount + 1;
  end
  clusterDistributions = clusterDistributions(validIdx);
  refinedClusterSamples = cat(1, refinedClusterSamples{:});
  refinedCliques = cat(1, refinedCliques{:});
  clusterThresholds = cat(1, clusterThresholds{:});
  refinedInstancePositions = cat(1, refinedInstancePositions{:});
  refinedClusters = cat(1, refinedClusters{:});
  eliminatedSamples = cat(1, eliminatedSamples{:});
  shownSamples = cat(1, refinedClusterSamples, eliminatedSamples);
  if ~isempty(shownSamples)
       shownClusters = cat(1, refinedClusters, ones(size(eliminatedSamples,1),1) * (max(refinedClusters) + 1));
  else
       shownClusters = [];
  end
end

function covMat = calculateCov(samples, dummyStd)
    sampleCov = eye(size(samples,2), 'single') * dummyStd;
    if size(samples,1) > 1
        covMat = cov(samples);
        covMat(1:(size(covMat,1) + 1):(size(covMat,1)^2)) = max(covMat(1:(size(covMat,1) + 1):(size(covMat,1)^2)), dummyStd);
    else
        covMat = sampleCov;
    end
    covMat = nearestSPD(covMat);
end