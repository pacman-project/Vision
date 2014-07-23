% the function computes element given clusterX aand clusterY

function [clusterNumber] = compute2elementIndex(clusterX, clusterY, nClusters)
    
    if clusterX <=0 || clusterY <=0 
        disp('ERRRROOORR!!');
    end
    clusterNumber = (clusterX - 1)* nClusters + clusterY;
    
end

