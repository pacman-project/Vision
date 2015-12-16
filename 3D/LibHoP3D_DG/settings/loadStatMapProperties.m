% This are the properties of the statistical map

function [statMapPropertiesAll, pairClusteringOptionsAll] = loadStatMapProperties(dataSetNumber)

    if dataSetNumber == 5  
        %% layer 3
        statMapProperties.vectStep = 5;
        statMapProperties.multX = 0.30;
        statMapProperties.multY = 0.20;
        statMapProperties.offsetMaximalAngle = 120;
        statMapProperties.angleStep = 0.01; % for saving file
        statMapProperties.xyzStep = 0.001;  % for saving file
        statMapPropertiesAll{3} = statMapProperties;
        
        pairClusteringOptions.penaltyThresh = 0.007;
        pairClusteringOptions.nCl_max = 11;
        pairClusteringOptions.alpha = 0.1;
        pairClusteringOptions.clusteringMethod = 3;
        pairClusteringOptions.sieveThresh = 100;
        pairClusteringOptionsAll{3} = pairClusteringOptions;
        
        %% layer 4
        statMapProperties.vectStep = 5;
        statMapProperties.multX = 0.30;
        statMapProperties.multY = 0.30;
        statMapProperties.offsetMaximalAngle = 100;
        statMapProperties.angleStep = 0.01; % for saving file
        statMapProperties.xyzStep = 0.001;  % for saving file
        statMapPropertiesAll{4} = statMapProperties;
        
        pairClusteringOptions.penaltyThresh = 0.007;
        pairClusteringOptions.nCl_max = 25;
        pairClusteringOptions.alpha = 0.1;
        pairClusteringOptions.clusteringMethod = 3;
        pairClusteringOptions.sieveThresh = 50;
        pairClusteringOptionsAll{4} = pairClusteringOptions;
        
        
        %% layer 5
        statMapProperties.vectStep = 5;
        statMapProperties.multX = 0.40;
        statMapProperties.multY = 0.10;
        statMapProperties.offsetMaximalAngle = 100;
        statMapProperties.angleStep = 0.01; % for saving file
        statMapProperties.xyzStep = 0.003;  % for saving file
        statMapPropertiesAll{5} = statMapProperties;
        
        pairClusteringOptions.penaltyThresh = 0.02;
        pairClusteringOptions.nCl_max = 25;
        pairClusteringOptions.alpha = 0.1;
        pairClusteringOptions.clusteringMethod = 3;
        pairClusteringOptions.sieveThresh = 25;
        pairClusteringOptionsAll{5} = pairClusteringOptions;
        
        %% layer 6
        statMapProperties.vectStep = 5;
        statMapProperties.multX = 0.40;
        statMapProperties.multY = 0.40;
        statMapProperties.offsetMaximalAngle = 100;
        statMapProperties.angleStep = 0.01; % for saving file
        statMapProperties.xyzStep = 0.003;  % for saving file
        statMapPropertiesAll{6} = statMapProperties;
        
        pairClusteringOptions.penaltyThresh = 0.02;
        pairClusteringOptions.nCl_max = 25;
        pairClusteringOptions.alpha = 0.1;
        pairClusteringOptions.clusteringMethod = 3;
        pairClusteringOptions.sieveThresh = 25;
        pairClusteringOptionsAll{6} = pairClusteringOptions;
    end

end

