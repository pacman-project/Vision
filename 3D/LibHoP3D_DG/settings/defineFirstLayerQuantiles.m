function [cluster1Centres, cluster1Bounds, thresh] = defineFirstLayerQuantiles(nClusters, dataSetNumber)
    
    if nClusters == 7
        if dataSetNumber > 0 && dataSetNumber <= 4 
            cluster1Bounds = [-40, -28, -15, -6, 6, 15, 28, 40];
            thresh = 40;
        elseif  dataSetNumber == 5
            cluster1Bounds = [-50, -35, -18, -8, 8, 18, 35, 50];
            thresh = 40;
        end
    end
    
    if nClusters == 9
        if dataSetNumber == 5
            cluster1Bounds = [-70, -50, -35, -20, -7, 7, 20, 35, 50, 70];
            thresh = 40;
        elseif dataSetNumber == 3
            cluster1Bounds = [-75   -50   -35   -20     -10    10    20    35    50    75];
            thresh = 40;
        end
    end
    
    if nClusters == 11
        if dataSetNumber == 5
            cluster1Bounds = [-90, -70, -50, -35, -20, -7, 7, 20, 35, 50, 70, 90];
            thresh = 40;
        elseif dataSetNumber == 3
            cluster1Bounds = [-78   -58   -39   -25   -15     -7    7    15    25    39    58    78];
            thresh = 40;
        end
    end
    
    
    %% compute cluster1Centres
    cluster1Centres = zeros(1, nClusters);
    
    for i = 1:nClusters
        cluster1Centres(i) = (cluster1Bounds(i) + cluster1Bounds(i+1))/2;
    end
    a = 2;

end

