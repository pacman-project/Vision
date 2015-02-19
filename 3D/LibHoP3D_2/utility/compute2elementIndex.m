% the function computes element given clusterX aand clusterY

function [clusterNumber] = compute2elementIndex(clusterX, clusterY, nClusters)
    
    if clusterX <=0 || clusterY <=0 
        disp('ERRRROOORR!!');
    end
    if clusterX == nClusters + 1 && clusterY == nClusters + 1  % this is an empty cell
        clusterNumber = nClusters^2 + 1;  % eg. 50 if nClusters == 7
        return;
    end
    
    
    clusterNumber = (clusterX - 1)* nClusters + clusterY;  
end

