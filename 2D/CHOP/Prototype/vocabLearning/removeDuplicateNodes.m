%> Name: removeDuplicateSubs
%>
%> Description: This function eliminates duplicate realizations from
%> graphLevel by checking whether the set of children of each part is unique.
%> Parts in vocabLevel which do not have any realizations after this
%> elimination are also deleted.
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
    graphNodeChildren = {graphLevel.children};
    sortedGraphNodeChildren = cellfun(@(x) mat2str(sort(x)), graphNodeChildren, 'UniformOutput', false);
    [~, validIdx, ~] = unique(sortedGraphNodeChildren, 'stable');
    graphLevel = graphLevel(validIdx);

    % After the realizations are removed, delete vocabulary nodes which
    % do not have any children left.
    allGraphLabelIds = [graphLevel.labelId];
    graphLabelIds = unique(allGraphLabelIds);
    bestSubs = bestSubs(graphLabelIds);
    vocabLevel = vocabLevel(graphLabelIds);
    for vocabNodeItr = 1:numel(vocabLevel)
        vocabLevel(vocabNodeItr).label = num2str(vocabNodeItr);
    end

    % Update graphLevel so that labelId of each node links to the
    % correct vocabulary part (now that some are deleted)
    graphLabelAssgnArr = zeros(max(graphLabelIds),1);
    graphLabelAssgnArr(graphLabelIds) = 1:numel(graphLabelIds);
    graphLabelAssgnArr = num2cell(graphLabelAssgnArr(allGraphLabelIds));
    [graphLevel.labelId] = deal(graphLabelAssgnArr{:});
end