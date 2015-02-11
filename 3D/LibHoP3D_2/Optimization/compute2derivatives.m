% the function computes clusterX and clusterY given number of the element

function [clusterX, clusterY] = compute2derivatives(elementNumber, nClusters)

    if mod(elementNumber, nClusters) == 0
        clusterY = nClusters;
        clusterX = round(elementNumber / nClusters);
    else
        clusterY = mod(elementNumber, nClusters);
        clusterX = ceil(elementNumber/nClusters);
    end
    
    if clusterX == nClusters + 1 && clusterY == 1 % this is an empty cell
        clusterY = nClusters + 1;     % ex. we return [8,8] if nClusters == 7
    end
    
end
