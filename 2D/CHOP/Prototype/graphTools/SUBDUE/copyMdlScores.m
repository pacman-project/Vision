function [childSubsFinal, childSubsExtend] = copyMdlScores(childSubsFinal, childSubsExtend)
    for subItr = 1:numel(childSubsFinal)
        childSubsExtend(subItr).mdlScore = childSubsFinal(subItr).mdlScore;
        childSubsExtend(subItr).normMdlScore = childSubsFinal(subItr).normMdlScore;
    end
end

