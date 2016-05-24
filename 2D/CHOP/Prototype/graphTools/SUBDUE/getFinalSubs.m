%> Name: getFinalSubs
%>
%> Description: Eliminates instances in each sub of childSubs that do not 
%> match within the limits specified by adaptiveThreshold.
%> 
%> @param childSubs List of subs to be processed.
%> @param adaptiveThreshold The threshold which sets a hard limit on
%> match costs of instances.
%>
%> @retval childSubsFinal Processed sub list. The subs are finalized.
%> @retval childSubs Pairs of childSubsFinal, with no adaptive threshold
%> limit applied. To be returned for further extensions.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 24.05.2014
function [childSubsFinal, childSubs] = getFinalSubs(childSubs, adaptiveThreshold)
    childSubsFinal = childSubs;
    validSubs = ones(numel(childSubs),1) > 0;
    parfor childSubItr = 1:numel(childSubs)
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
        childSubsFinal(childSubItr).instanceSigns = childSubsFinal(childSubItr).instanceSigns(validInstanceIdx, :);
        childSubsFinal(childSubItr).instanceValidationIdx = childSubsFinal(childSubItr).instanceValidationIdx(validInstanceIdx, :);
    end
    childSubsFinal = childSubsFinal(validSubs);
    childSubs = childSubs(validSubs);
end