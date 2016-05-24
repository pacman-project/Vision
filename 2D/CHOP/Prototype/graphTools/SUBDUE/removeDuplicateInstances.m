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
        [~, validIdx, ~] = unique(sort(instanceChildren,2), 'rows', 'stable');

        % Get minimum matching costs and children.
        instanceChildren = instanceChildren(validIdx, :);

        % Finally, order children by rows and update remaining data
        % structures.
        [subList(bestSubItr).instanceChildren, idx] = sortrows(instanceChildren);
        validIdx = validIdx(idx);
        subList(bestSubItr).instanceCenterIdx = subList(bestSubItr).instanceCenterIdx(validIdx, :);
        if ~isempty(subList(bestSubItr).instanceEdges)
            subList(bestSubItr).instanceEdges = subList(bestSubItr).instanceEdges(validIdx, :);
        end
        subList(bestSubItr).instanceCategories = subList(bestSubItr).instanceCategories(validIdx, :);
        subList(bestSubItr).instanceSigns = subList(bestSubItr).instanceSigns(validIdx, :);
        subList(bestSubItr).instanceValidationIdx = subList(bestSubItr).instanceValidationIdx(validIdx, :);
        subList(bestSubItr).instanceMappings = subList(bestSubItr).instanceMappings(validIdx, :);
    end
end

