%> Name: removeDuplicateInstances
%>
%> Description: This function removes duplicate instances from each sub in
%> subList, and orders the instances based on their match costs.
%> 
%> @param subList List of subs to be processed.
%>
%> @retval subList The list of subs with eliminated instances.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 23.09.2015
function [ subList ] = removeDuplicateInstances( subList )
    parfor bestSubItr = 1:numel(subList)
        instanceChildren = subList(bestSubItr).instanceChildren;
        instanceMatchCosts = subList(bestSubItr).instanceMatchCosts;
        % We handle these cases by only keeping
        % unique instances. In addition, for each instance, the minimum
        % cost of matching is kept here.
        [minMatchCosts, sortIdx] = sort(instanceMatchCosts, 'ascend');
        sortedChildren = instanceChildren(sortIdx, :);
        [~, validIdx, ~] = unique(sort(sortedChildren,2), 'rows', 'stable');
        sortIdx = sortIdx(validIdx);

        % Get minimum matching costs and children.
        instanceMatchCosts = minMatchCosts(validIdx, :);
        sortedChildren = sortedChildren(validIdx, :);

        % Finally, order children by rows and update remaining data
        % structures.
        [subList(bestSubItr).instanceChildren, idx] = sortrows(sortedChildren);
        sortIdx = sortIdx(idx);
        subList(bestSubItr).instanceMatchCosts = instanceMatchCosts(idx, :);
        subList(bestSubItr).instanceCenterIdx = subList(bestSubItr).instanceCenterIdx(sortIdx, :);
        if ~isempty(subList(bestSubItr).instanceEdges)
            subList(bestSubItr).instanceEdges = subList(bestSubItr).instanceEdges(sortIdx, :);
        end
        subList(bestSubItr).instanceCategories = subList(bestSubItr).instanceCategories(sortIdx, :);
        subList(bestSubItr).instanceSigns = subList(bestSubItr).instanceSigns(sortIdx, :);
        subList(bestSubItr).instanceValidationIdx = subList(bestSubItr).instanceValidationIdx(sortIdx, :);
        subList(bestSubItr).instanceMappings = subList(bestSubItr).instanceMappings(sortIdx, :);
        subList(bestSubItr).instanceExactMatchFlags = subList(bestSubItr).instanceExactMatchFlags(sortIdx, :);
    end
end

