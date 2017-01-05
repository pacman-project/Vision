function [partClusters, refinedClusterSamples, refinedCliques, refinedClusters, refinedInstancePositions, shownSamples, shownClusters, clusterDistributions, numberOfClusters] = learnJointStats(jointPositions, instancePositions, savedCliques, realRFSize, maxClusters) 
  % We will experiment with different algorithms here, 
  % from completely unsupervised to supervised.
  dummyStd = 0.001;
  stdThr = 2.576;
  maxPointsToCluster = 100;
  if size(jointPositions,1) == 1
       partClusters = 1;
  else
       % Select a number of rows and cluster them, not all data.
       selectedRows = datasample(1:size(jointPositions,1), min(size(jointPositions,1), maxPointsToCluster), 'Replace', false);
       selectedPositions = jointPositions(selectedRows, :);
       [selectedClusters, ~] = gmeans(double(selectedPositions*(1/((realRFSize-1) * sqrt(size(selectedPositions,2))))), maxClusters, 0.05, @checkGaussian);
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
  refinedClusterSamples = cell(numberOfClusters,1);
  refinedInstancePositions = cell(numberOfClusters,1);
  refinedClusters = cell(numberOfClusters,1);
  refinedCliques = cell(numberOfClusters,1);
  eliminatedSamples = cell(numberOfClusters,1);
  for clusterItr = 1:numberOfClusters
      idx = partClusters == clusterItr;
      clusterSamples = jointPositions(idx,:);
      clusterPositions = instancePositions(idx, :);
      clusterCliques = savedCliques(idx,:);
      clusterCenters(clusterItr,:) = sum(clusterSamples,1) / nnz(idx);

      % We eliminate nodes that are far away from the cluster centers
      % in every cluster.
      covMat = calculateCov(clusterSamples, dummyStd);
      distances = pdist2(clusterCenters(clusterItr,:), clusterSamples, 'mahalanobis', covMat);

      % Filter out unlikely instances and save the data.
      thresh = max(min(distances)+0.0001, stdThr); % 98 pct confidence interval
      validInstances = distances <= thresh;
      refinedClusterSamples{clusterItr} = clusterSamples(validInstances,:);
      refinedCliques{clusterItr} = clusterCliques(validInstances, :);
      refinedInstancePositions{clusterItr} = clusterPositions(validInstances,:);
      newCovMat = calculateCov(clusterSamples(validInstances,:), dummyStd);
      eliminatedSamples{clusterItr} = clusterSamples(~validInstances, :);
      refinedClusters{clusterItr} = ones(nnz(validInstances),1) * clusterItr;

      % Save covariance matrix.
      clusterDistributions{clusterItr} = gmdistribution(clusterCenters(clusterItr,:), newCovMat);
  end
  refinedClusterSamples = cat(1, refinedClusterSamples{:});
  refinedCliques = cat(1, refinedCliques{:});
  refinedInstancePositions = cat(1, refinedInstancePositions{:});
  refinedClusters = cat(1, refinedClusters{:});
  eliminatedSamples = cat(1, eliminatedSamples{:});
  shownSamples = cat(1, refinedClusterSamples, eliminatedSamples);
  shownClusters = cat(1, refinedClusters, ones(size(eliminatedSamples,1),1) * (max(refinedClusters) + 1));
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