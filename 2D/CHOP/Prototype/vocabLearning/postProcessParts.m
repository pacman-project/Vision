function [vocabLevel, redundantVocabLevel, graphLevel, newDistanceMatrix] = postProcessParts(vocabLevel, subClasses, newDistanceMatrix, graphLevel)
    % Assign new labels of the remaining realizations.
    [remainingComps, ~, IC] = unique([graphLevel.labelId]);
    IC = num2cell(int32(IC));
    [graphLevel.labelId] = deal(IC{:});

    % Get the redundant compositions to help inference process.
    if ~isempty(subClasses)
        redundantComps = ismember(subClasses, remainingComps) & subClasses ~= (1:numel(subClasses))';
        redundantVocabLevel = vocabLevel(1, redundantComps);
    else
        redundantVocabLevel = [];
    end
    
    % Assign correct labels to the redundant compositions.
    if ~isempty(redundantVocabLevel)
        maxRemaining = max(remainingComps);
        remainingMatchArr = zeros(maxRemaining,1, 'int32');
        remainingMatchArr(remainingComps) = 1:numel(remainingComps);
        VIC = num2cell(remainingMatchArr(subClasses(redundantComps)));
        [redundantVocabLevel.label] = deal(VIC{:});
    end
    
    % Eliminate unused compositions from vocabulary.
    vocabLevel = vocabLevel(1, remainingComps);
    newLabelArr = num2cell(int32(1:numel(vocabLevel)));
    [vocabLevel.label] = deal(newLabelArr{:});
    
    % Assign distance matrix.
    if ~isempty(newDistanceMatrix)
        numberOfRemainingComps = numel(remainingComps);
        newDistanceMatrix = newDistanceMatrix(remainingComps, remainingComps);

        % Normalize distances by the size of compared parts.
        childrenCounts = {vocabLevel.children};
        childrenCounts = cellfun(@(x) numel(x), childrenCounts);
        for partItr = 1:numel(vocabLevel);
            newDistanceMatrix(partItr,:) = newDistanceMatrix(partItr,:) ./ ...
               (max(childrenCounts, repmat(childrenCounts(partItr), 1, numberOfRemainingComps)) * 2 - 1);
        end
    end
    newDistanceMatrix = newDistanceMatrix / max(max(newDistanceMatrix));
end

