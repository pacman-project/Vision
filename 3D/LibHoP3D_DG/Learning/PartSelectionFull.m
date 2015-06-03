% this is a script to perform full part selection procedure (iterative procedure, greedy algorithm)

function [triplesCurOut, coverageOut, nCurClusters] = PartSelectionFull( nClusters, n2Clusters, dataSetNumber, fieldSize, iterations,...
                             layerID, fileForVisualizationPrevLayer, lenSelected, cluster1Centres, numSimilar, offsetsConventional, depthStep, ...
                             is_multyScale, scales, lineAdders, pathBase, is_subset, subsetPercent, root, dsN, nCl)
                         
    for i = length(scales):-1:1 % going from the smallest scale to the largest one!
        
        input_path = getPathScale(pathBase, lineAdders{i});  
        [list_input, ~, ~, lenF] = extractFileListGeneral(input_path, is_subset, subsetPercent, dataSetNumber);  
        checkImages(list_input, lenF);   % makes images in the folder 3 channels ones
        
        
        [~, ~, ~, statisticsLayerSieved, statisticsLayerAggregated] = GetStatisticsFiles(root, dsN, nCl, lineAdders{i}, layerID);
         
        [triplesCurOut, coverageOut, ~] = partSelectionNew(nClusters, n2Clusters, statisticsLayerSieved, statisticsLayerAggregated, ...
             dataSetNumber, fieldSize, list_input, lenF, iterations, layerID, fileForVisualizationPrevLayer, cluster1Centres, ...
             lenSelected, numSimilar, offsetsConventional, depthStep);
         
         
         if i > 1 % project the learned parts to the higher scale in order not to repeat learning
             reprojectCoverage(list_input, lineAdders{i}, lineAdders{i-1});
         end       
        
    end



    [coverageOut, inds] = sort(coverageOut, 'descend');
    triplesCurOut = triplesCurOut(inds, :);

    % select only parts with large enough coverage
    idx = coverageOut(coverageOut > lenSelected);
    nCurClusters = length(idx);
    
end

