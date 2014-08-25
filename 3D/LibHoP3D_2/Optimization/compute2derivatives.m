% the function computes clusterX and clusterY given number of the element

function [clusterX, clusterY] = compute2derivatives(elementNumber, nClusters)

    if mod(elementNumber, nClusters) == 0
        clusterY = nClusters;
        clusterX = round(elementNumber / nClusters);
    else
        clusterY = mod(elementNumber, nClusters);
        clusterX = ceil(elementNumber/nClusters);
    end
    
end
