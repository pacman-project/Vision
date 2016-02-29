% This are the properties of the statistical map

function [statMapPropertiesAll, pairClusteringOptionsAll] = loadStatMapProperties(dataSetNumber)

    if dataSetNumber == 5  
        % layer:     3     4     5     6      7      8      9     10
        quaternionSteps = {0.02, 0.02, 0.04,  0.04,  0.06,  0.06, 0.08,  0.08};  % for statistical maps
        multXs   = {.22,  .2,   0.15, 0.13,   0.15,  0.15, 0.11,  0.15};
        multYs   = {.15,  .2,   0.1,  0.13,   0.15   0.15, 0.13,  0.13};
%         offMA    = {120,  100,  100,  100,    100,   100};
        angleSteps={.01,  0.01, 0.01, 0.01,   0.01, 0.01,  0.01,  0.01};                  % for file storage
    xyzSteps = {0.00066, 0.00066, 0.003, 0.003, 0.01, 0.01, 0.03, 0.03};
    penTreshes={0.016,   0.04,    0.05, 0.1,   0.1,  0.15,  0.1,  0.15};
        sieveThs = {20, 100,  100,    100,    100,    50,   50,   50};
%         alphas   = {0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12};
    alphas = {     0.8,   0.4,  0.4,   0.4,  0.6,  0.6,   0.6,  0.5,  1,1,1};
        b2wRatio = {3,    20,     3,     3,    3,    3,   3,  3,   3};
        
        nCl_max = {13,   13,     13,    13,    13,    25,  25,  25};
        clusteringMethod = 3;
        
        
%         
%         %% layer 7
%         statMapProperties.vectStep = 10;
%         statMapProperties.multX = 0.2; % becomes multD
%         statMapProperties.multY = 0.10;
%         statMapProperties.offsetMaximalAngle = 100;
%         statMapProperties.angleStep = 0.02; % for saving file
%         statMapProperties.xyzStep = 0.01;  % for saving file
%         statMapPropertiesAll{7} = statMapProperties;
%         
%         pairClusteringOptions.penaltyThresh = 0.1;
%         pairClusteringOptions.nCl_max = 9;
%         pairClusteringOptions.alpha = 0.1;
%         pairClusteringOptions.clusteringMethod = 3;
%         pairClusteringOptions.sieveThresh = 10;
%         pairClusteringOptions.b2wRatio = 3;
%         pairClusteringOptionsAll{7} = pairClusteringOptions;
    end
    
    for layer = 3:8
            
        i = layer - 2;
        statMapProperties.quaternionSteps = quaternionSteps{i};
        statMapProperties.multX = multXs{i};
        statMapProperties.multY = multYs{i};
        statMapProperties.offsetMaximalAngle = 1.0;   % maximal value of quaternions         %offMA{i};
        statMapProperties.angleStep = angleSteps{i}; % for saving file
        statMapProperties.xyzStep = xyzSteps{i};  % for saving file
        statMapPropertiesAll{layer} = statMapProperties;

        pairClusteringOptions.penaltyThresh = penTreshes{i};
        pairClusteringOptions.nCl_max = nCl_max{i};
        pairClusteringOptions.alpha = alphas{i};
        pairClusteringOptions.clusteringMethod = clusteringMethod;
        pairClusteringOptions.sieveThresh = sieveThs{i};
        pairClusteringOptions.b2wRatio = b2wRatio{i};
        pairClusteringOptionsAll{layer} = pairClusteringOptions;
    end

end

