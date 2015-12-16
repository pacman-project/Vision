% converts a part to surface descriptor

function [X_first, nPrevClusters, emptyIndicator] = convertToSurfaceDescriptor(X, lenCombs, layerID, nClusters, n2Clusters, fileForVisualizationPrevLayer,  ...
                       offsetsConventional, depthStep, fieldCenter, cluster1Centres, downsamplingScheme, clusterCurDepths)

    disp('computing surface descriptors');
                   
    for i = 3:8
        tripleOutDepth{i} = [];
    end               
                   
    if layerID > 3
        load(fileForVisualizationPrevLayer);
    end
      
    if layerID == 2
        nPrevClusters = nClusters;
        tripleOutDepth = [];
    elseif layerID == 3
        nPrevClusters = n2Clusters;
        tripleOutDepth{3} = store3Layer(X, clusterCurDepths, lenCombs, nClusters);
    else
       nPrevClusters = size(tripleOutDepth{layerID - 1}, 1);
       tripleOutDepth{layerID} = store4Layer(X, clusterCurDepths, lenCombs, nClusters);
    end
        
    
    X_first = [];
    emptyIndicator = [];
    
    parfor i = 1:lenCombs    
        
        [positions, elements] = partMeanReconstruction(layerID, i, fieldCenter, tripleOutDepth, offsetsConventional, depthStep, nClusters);
                                                                                                           
        [X_first_lile, emptyIndicatorLine] = partToSurfaceDescriptor(elements, positions, nClusters, cluster1Centres, downsamplingScheme);
        X_first = [X_first; X_first_lile];
        emptyIndicator = [emptyIndicator; emptyIndicatorLine];
    end

end

