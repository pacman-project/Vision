%> Name: discoverSubs
%>
%> Description: This functions discovers hidden structures in graphLevel,
%> and creates nodes of the new graph that is a compressed version of
%> graphLevel. Please note that the number of nodes in output graphLevel
%> could be higher than the input graphLevel, which does not make much sense
%> from a compression standpoint. However, there is a lot of redundancy in
%> the output. This problem is dealt with later in the inhibition step.
%>
%> @param vocabLevel If preDefinedSearch is 1, the compositions in this
%> vocabulary level are searched in graphLevel. If 0, simply ignored.
%> @param graphLevel The current object graphs' level.
%> @param options Program options.
%> @param currentFolder Path to the workspace folder.
%> @param preDefinedSearch If true, supervised search is run. 
%> If empty, unsupervised SUBDUE runs over graphLevel.
%> @param levelItr current level id.
%>
%> @retval vocabLevel If preDefinedSearch is 1, [] is returned. If it is 0,
%> best compositions discovered are returned.
%> @retval graphLevel The graph consisting of part realizations discovered.
%> It's just a dummy graph including realization labels and children.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 15.01.2014
%> 'self' type search added on 05.02.2014
function [vocabLevel, graphLevel, optimalThreshold, isSupervisedSelectionRunning, previousAccuracy] = discoverSubs( vocabLevel, graphLevel, ...
    nodeDistanceMatrix, options, threshold, levelItr, supervisedSelectionFlag, isSupervisedSelectionRunning, previousAccuracy, level1Nodes)
    startTime = tic;
    edgeDistanceMatrix = options.edgeDistanceMatrix;
    display(['.... Discovering compositions in level ' num2str(levelItr) '.']); 
    load([options.currentFolder '/output/' options.datasetName '/export.mat'], 'categoryArrIdx', 'validationIdx');
    
    % Search for substructures.
    [vocabLevel, graphLevel, optimalThreshold, isSupervisedSelectionRunning, previousAccuracy] = runSubdue(vocabLevel, graphLevel, threshold,...
        nodeDistanceMatrix, edgeDistanceMatrix, categoryArrIdx, validationIdx, ...
        supervisedSelectionFlag, isSupervisedSelectionRunning, previousAccuracy, level1Nodes, options);
    
    % Show time elapsed.
    display(['.... Time elapsed: ' num2str(toc(startTime)) ' secs.']);
    display(['.... Found ' ...
        num2str(numel(graphLevel)) ' instances of ' num2str(numel(vocabLevel)) ' compositions.']);
end
