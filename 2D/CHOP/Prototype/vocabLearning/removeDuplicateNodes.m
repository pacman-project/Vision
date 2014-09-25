%> Name: removeDuplicateSubs
%>
%> Description: This function eliminates duplicate realizations from
%> graphLevel by checking whether the graph descriptions of two subs are
%> essentially the same. If so, only one such part is kept in the
%> vocabulary.
%>
%> @param vocabLevel Input set of parts.
%> @param graphLevel Input set of realizations.
%> @param bestSubs List of parts in terms of substructures.
%>
%> @retval vocabLevel Output (valid) set of parts.
%> @retval graphLevel Unique (in terms of children) set of realizations.
%> @param bestSubs Unique list of parts in terms of substructures.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 19.09.2014
function [vocabLevel, graphLevel, bestSubs] = removeDuplicateNodes(vocabLevel, graphLevel, bestSubs)
    % Sanity check on params.
    if isempty(graphLevel) || isempty(vocabLevel) || isempty(bestSubs)
       return; 
    end
    
    % Re-order edges of bestSubs in order to get a unique description that
    % is invariant of the inner ordering of edges.
    subEdges = cell(1, numel(bestSubs));
    for subItr = 1:numel(bestSubs)
        subEdges(subItr) = {sortrows(bestSubs(subItr).edges)};
    end
    subCenters = {bestSubs.centerId};
    vocabDescriptors = cellfun(@(x,y) num2str([x; y(:)]'), subCenters, subEdges, 'UniformOutput', false);
    [~, validVocabLevelIdx, ~] = unique(vocabDescriptors, 'stable');
    
    % Eliminate duplicate parts and their realizations.
    vocabLevel = vocabLevel(validVocabLevelIdx);
    validGraphLevelIdx = ismember(cat(1,graphLevel.labelId), validVocabLevelIdx); 
    graphLevel = graphLevel(validGraphLevelIdx);

    % After the realizations are removed, delete vocabulary nodes which
    % do not have any children left.
    allGraphLabelIds = [graphLevel.labelId];
    bestSubs = bestSubs(validVocabLevelIdx);
    for vocabNodeItr = 1:numel(vocabLevel)
        vocabLevel(vocabNodeItr).label = num2str(vocabNodeItr);
    end

    % Update graphLevel so that labelId of each node links to the
    % correct vocabulary part (now that some are deleted)
    graphLabelAssgnArr = zeros(max(validVocabLevelIdx),1);
    graphLabelAssgnArr(validVocabLevelIdx) = 1:numel(validVocabLevelIdx);
    graphLabelAssgnArr = num2cell(graphLabelAssgnArr(allGraphLabelIds));
    [graphLevel.labelId] = deal(graphLabelAssgnArr{:});
end