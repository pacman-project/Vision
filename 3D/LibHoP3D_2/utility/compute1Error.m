

function [error] = compute1Error(clusterX, cluster1Centres, cluster1Lengths, nClusters)
    
    if clusterX <= cluster1Centres(1)
        error = abs(  2 *(clusterX - cluster1Centres(1))/cluster1Lengths(1)  );
    elseif clusterX >= cluster1Centres(nClusters)
        error = abs(  2 * (clusterX - cluster1Centres(nClusters))/cluster1Lengths(nClusters)  );
    else 
        error = abs(  (clusterX - cluster1Centres(clusterI)) / cluster1Lengths(clusterI)  );
    end
    
    if error > 1
        error = 1;
    end

end

