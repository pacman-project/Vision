function pairClusteringOptionsAll = loadpairClusteringOptions(dataSetNumber, offsetConventional)

    if dataSetNumber == 5
        
        % layerID == 3;
        pairClusteringOptions.penaltyThresh = offsetConventional{3};
        pairClusteringOptions.nCl_max = 45;
        pairClusteringOptions.alpha = 0.11;
        pairClusteringOptions.clusteringMethod = 3; % hierarchical clustering
        pairClusteringOptions.sieveThresh = 100;
        pairClusteringOptionsAll{3} = pairClusteringOptions;
        
        % layerID == 4
        
    end


end

