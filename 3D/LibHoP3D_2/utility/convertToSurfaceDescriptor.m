% converts a part to surface descriptor

function [X_first, nPrevClusters, emptyIndicator] = convertToSurfaceDescriptor(X, lenCombs, layerID, nClusters, n2Clusters, fileForVisualizationPrevLayer,  ...
                       displ3, displ5, displ7, fieldCenter, cluster1Centres, downsamplingScheme, clusterCurDepths)

    for i = 3:8
        tripleOutDepth{i} = [];
    end               
                   
    if layerID > 3
        load(fileForVisualizationPrevLayer);
    end
        
    if layerID == 3
        nPrevClusters = n2Clusters;
        tripleOutDepth{3} = store3Layer(X, clusterCurDepths, lenCombs, nClusters, 1);
    else
       nPrevClusters = size(tripleOutDepth{layerID - 1}, 1);
       tripleOutDepth{layerID} = store4Layer(X, clusterCurDepths, lenCombs, nClusters, 1);
    end
        
    
    X_first = [];
    emptyIndicator = [];
    
    parfor i = 1:lenCombs    
        
        [positions, elements] = partMeanReconstruction(layerID, i, fieldCenter, tripleOutDepth, displ3, displ5, displ7, nClusters);
                                                                                                           
        [X_first_lile, emptyIndicatorLine] = partToSurfaceDescriptor(elements, positions, nClusters, cluster1Centres, downsamplingScheme);
        X_first = [X_first; X_first_lile];
        emptyIndicator = [emptyIndicator; emptyIndicatorLine];
    end

end

