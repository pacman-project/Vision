% define centres and sizes of the first layer parts

function [clusterCenters, outClusterSize] = defineCluster1Centers(nClusters, clusterSizes, thresh)

     outClusterSize = clusterSizes * thresh * 2;     
     clusterCenters = zeros(1, nClusters);
    
     curShift = 0;
     
     for i = 1:nClusters
        clusterCenters(i) = -thresh + curShift + outClusterSize(i)/2;
        curShift = curShift + outClusterSize(i);
     end
     
end

