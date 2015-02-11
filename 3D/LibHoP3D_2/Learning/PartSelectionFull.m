% this is a script to perform full part selection procedure (iterative procedure, greedy algorithm)

function [triplesCurOut, coverageOut, nCurClusters] = PartSelectionFull( nClusters, n2Clusters, statisticsPrevLayerSieved, statisticsPrevLayerAggregated, ...
                             dataSetNumber, fieldSize, list_depth, lenF, iterations,...
                             layerID, fileForVisualizationPrevLayer, lenSelected, displ3, displ5, displ7, cluster1Centres, depthStep, numSimilar, is_GPU_USED)

    checkImages(list_depth, lenF);   % makes images in the folder 3 channels ones

    [triplesCurOut, coverageOut, ~] = partSelectionNew(nClusters, n2Clusters, statisticsPrevLayerSieved, statisticsPrevLayerAggregated, ...
                             dataSetNumber, fieldSize, list_depth, lenF, ...
                             iterations, layerID, fileForVisualizationPrevLayer, displ3, displ5, displ7, cluster1Centres, ...
                             depthStep, lenSelected, numSimilar, is_GPU_USED);

    [coverageOut, inds] = sort(coverageOut, 'descend');
    triplesCurOut = triplesCurOut(inds, :);

    % select only parts with large enough coverage
    idx = coverageOut(coverageOut > lenSelected);
    nCurClusters = length(idx);
    
end

