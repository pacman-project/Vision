function [ vocabLevel, graphLevel ] = reorderSubs( vocabLevel, graphLevel )
    %% To start with, we reorder the graph level based on repetition of words.
    % This is done to alleviate effects of non-optimal mdl score
    % calculation of SUBDUE.
    if ~(numel(graphLevel) > 0 && numel(vocabLevel) > 0)
       graphLevel = [];
       vocabLevel = [];
       return; 
    end
    
    newGraphLevel = graphLevel;
    labelIds = [graphLevel.labelId];
    labelOccurences = zeros(max(labelIds),1);
    for labelItr = 1:max(labelIds)
       labelOccurences(labelItr) = numel(find(labelIds==labelItr)) * numel(vocabLevel(labelItr).children);
    end
    [~, idx] = sort(labelOccurences, 'descend');
    vocabLevel = vocabLevel(idx);
    
    % Reassign graph level by ordering them in based on frequency of words.
    realizationOffset = 1;
    for labelItr = 1:numel(idx)
        assignedRealizationIdx = find(labelIds==idx(labelItr));
        assignedRealizations = graphLevel(assignedRealizationIdx);
        assignedRealizations = arrayfun(@(s) setfield(s, 'labelId', labelItr),assignedRealizations);
        newGraphLevel(realizationOffset:(realizationOffset+numel(assignedRealizations)-1)) = ...
            assignedRealizations;
        realizationOffset = realizationOffset + numel(assignedRealizationIdx);
    end
    graphLevel = newGraphLevel;
end

