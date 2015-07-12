function [childSubsFinal, childSubs] = getFinalSubs(childSubs, adaptiveThreshold)
    childSubsFinal = childSubs;
    validSubs = ones(numel(childSubs),1) > 0;
    for childSubItr = 1:numel(childSubs)
        validInstanceIdx = childSubs(childSubItr).instanceMatchCosts < adaptiveThreshold;
        if nnz(validInstanceIdx) == 0
            validSubs(childSubItr) = 0;
            continue;
        end
        childSubsFinal(childSubItr).instanceChildren = childSubsFinal(childSubItr).instanceChildren(validInstanceIdx, :);
        childSubsFinal(childSubItr).instanceMappings = childSubsFinal(childSubItr).instanceMappings(validInstanceIdx, :);
        childSubsFinal(childSubItr).instanceCenterIdx = childSubsFinal(childSubItr).instanceCenterIdx(validInstanceIdx, :);
        childSubsFinal(childSubItr).instanceCategories = childSubsFinal(childSubItr).instanceCategories(validInstanceIdx, :);
        childSubsFinal(childSubItr).instanceEdges = childSubsFinal(childSubItr).instanceEdges(validInstanceIdx, :);
        childSubsFinal(childSubItr).instanceMatchCosts = childSubsFinal(childSubItr).instanceMatchCosts(validInstanceIdx, :);
        childSubsFinal(childSubItr).instanceSigns = childSubsFinal(childSubItr).instanceSigns(validInstanceIdx, :);
        childSubsFinal(childSubItr).instanceValidationIdx = childSubsFinal(childSubItr).instanceValidationIdx(validInstanceIdx, :);
    end
    childSubsFinal = childSubsFinal(validSubs);
    childSubs = childSubs(validSubs);
end