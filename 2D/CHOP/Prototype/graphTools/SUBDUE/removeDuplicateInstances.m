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
    subCenters = cat(1,subList.centerId);
    subEdges = {subList.edges};
    
    parfor bestSubItr = 1:numel(subList)
        instanceChildren = subList(bestSubItr).instanceChildren;
        if size(instanceChildren,2) > 2
             curEdges = subEdges{bestSubItr};
             if ismember(subCenters(bestSubItr), curEdges(:,2))
                  [~, validIdx, ~] = unique(sort(instanceChildren,2), 'rows', 'stable');
             else
                  validIdx = (1:size(instanceChildren,1))';
             end
        else
             validIdx = (1:size(instanceChildren,1))';
        end

        % Get minimum matching costs and children.
        instanceChildren = instanceChildren(validIdx, :);

        % Finally, order children by rows and update remaining data
        % structures.
        [subList(bestSubItr).instanceChildren, idx] = sortrows(instanceChildren);
        validIdx = validIdx(idx);
        subList(bestSubItr).instanceCenterIdx = subList(bestSubItr).instanceCenterIdx(validIdx, :);
        subList(bestSubItr).instanceSigns = subList(bestSubItr).instanceSigns(validIdx, :);
    end
end

