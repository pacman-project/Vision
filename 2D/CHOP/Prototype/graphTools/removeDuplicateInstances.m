function [ bestSubs ] = removeDuplicateInstances( bestSubs )
    parfor bestSubItr = 1:numel(bestSubs)
        instanceChildren = bestSubs(bestSubItr).instanceChildren;
        instanceMatchCosts = bestSubs(bestSubItr).instanceMatchCosts;
        % We handle these cases by only keeping
        % unique instances. In addition, for each instance, the minimum
        % cost of matching is kept here.
        [minMatchCosts, sortIdx] = sort(instanceMatchCosts, 'ascend');
        sortedChildren = instanceChildren(sortIdx, :);
        [~, validIdx, ~] = unique(sortedChildren, 'rows', 'stable');
        sortIdx = sortIdx(validIdx);

        % Get minimum matching costs and children.
        instanceMatchCosts = minMatchCosts(validIdx, :);
        sortedChildren = sortedChildren(validIdx, :);

        % Finally, order children by rows.
        [bestSubs(bestSubItr).instanceChildren, idx] = sortrows(sortedChildren);
        sortIdx = sortIdx(idx);
        bestSubs(bestSubItr).instanceMatchCosts = instanceMatchCosts(idx, :);
        bestSubs(bestSubItr).instanceCenterIdx = bestSubs(bestSubItr).instanceCenterIdx(sortIdx, :);
        if ~isempty(bestSubs(bestSubItr).instanceEdges)
            bestSubs(bestSubItr).instanceEdges = bestSubs(bestSubItr).instanceEdges(sortIdx, :);
        end
        bestSubs(bestSubItr).instanceCategories = bestSubs(bestSubItr).instanceCategories(sortIdx, :);
        bestSubs(bestSubItr).instanceSigns = bestSubs(bestSubItr).instanceSigns(sortIdx, :);
    end
end

