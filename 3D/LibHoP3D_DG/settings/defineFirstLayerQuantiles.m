function [cluster1Centres, cluster1Bounds, thresh] = defineFirstLayerQuantiles(nClusters, dataSetNumber)
    
    if nClusters == 7
        if dataSetNumber > 0 && dataSetNumber <= 4 
            cluster1Centres = [-34, -21.5, -10.5, 0, 10.5, 21.5, 34];
            cluster1Bounds = [-40, -28, -15, -6, 6, 15, 28, 40];
            thresh = 40;
        end
    end
    

end

