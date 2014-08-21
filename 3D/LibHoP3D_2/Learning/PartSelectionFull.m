% this is a script to perform full part selection procedure (iterative procedure, greedy algorithm)

function [triples3Out, coverageOut, n3Clusters] = PartSelectionFull( nClusters, n2Clusters, statisticsPrevLayerSieved, statisticsPrevLayerAggregated, ...
                             dataSetNumber, fieldSize, list_depth, lenF, meargeThresh, iterations,...
                             layerID, fileForVisualizationPrevLayer, lenSelected)

    checkImages(list_depth, lenF);   % makes images in the folder 3 channels ones

    [triples3Out, coverageOut, ~] = partSelectionNew(nClusters, n2Clusters, statisticsPrevLayerSieved, statisticsPrevLayerAggregated, ...
                             dataSetNumber, fieldSize, list_depth, lenF, meargeThresh, ...
                             iterations, layerID, fileForVisualizationPrevLayer);

    [coverageOut, inds] = sort(coverageOut, 'descend');
    triples3Out = triples3Out(inds, :);

    % select only parts with large enough coverage
    idx = coverageOut(coverageOut > lenSelected);
    n3Clusters = length(idx);
    
end

