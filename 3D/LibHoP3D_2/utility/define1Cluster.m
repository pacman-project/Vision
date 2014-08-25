% This function defines a cluster dx should belong, having parameters
% dx, thresh and nClusters

% now we pretend that nClusters = 9


% -thresh |------------------------------------------------------| thresh

function clusterX = define1Cluster(dX, cluster1Bounds, nClusters)

    clusterX = imquantize(dX, cluster1Bounds) - 1;
    clusterX(clusterX == nClusters + 1) = 0;
     
end

