% the function computes element given clusterX aand clusterY

function [clusterNumber] = compute2elementIndex(clusterX, clusterY, nClusters)

    clusterNumber = (clusterX - 1)* nClusters + clusterY;
    
end

