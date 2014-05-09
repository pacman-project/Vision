% This function defines a cluster dx should belong, having parameters
% dx, thresh and nClusters

% now we pretend that nClusters = 9


% -thresh |------------------------------------------------------| thresh

function clusterX = define1Cluster(dX, nClusters, clusterSizes, thresh)

     if dX <= -thresh +  clusterSizes(1)/2
        clusterX = 1;
     elseif dX >= thresh - clusterSizes(end)/2
        clusterX = nClusters;  
     else
          curCluster = 1;
          right = - thresh + (clusterSizes(1)/2);
          while dX >= right
              curCluster = curCluster + 1;
              right = right + (clusterSizes(curCluster));
          end
          clusterX = curCluster;
     end
end

