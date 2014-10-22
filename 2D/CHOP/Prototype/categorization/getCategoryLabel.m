%> Name: getCategoryLabel
%>
%> Description: Estimates a category label for the given realization set.
%> The method checks for highest level realizations, and tries to determine
%> whether they lean toward a certain category. In case it cannot determine,
%> it checks the lower level to get rid of ambiguities (recursively).
%>
%> @param vocabulary Learned vocabulary including category probabilities for 
%> each node.
%> @param exportArr Set of realizations.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 04.07.2014
function [ categoryLabel ] = getCategoryLabel(vocabulary, exportArr)
%    maxLevels = max(exportArr(:,4));
    maxLevels = 2;
    categoryLabel = -1;
    numberOfCategories = numel(vocabulary{1}(1).categoryArr);
    categoryList = 1:numberOfCategories;
    categoryDecisionArr = zeros(1, numberOfCategories, 'single');
    bestCategories = [];
    for levelItr = maxLevels:-1:1
        vocabLevel = vocabulary{levelItr};
        categoryArrs = cat(1, vocabLevel.categoryArr);
        if isempty(categoryArrs) 
           break; 
        end
        nodes = exportArr(exportArr(:,4) == levelItr,:);
        if size(nodes,1)>0
            categoryDecisionArr = categoryDecisionArr + sum(categoryArrs(nodes(:,1),:),1) / size(nodes,1);
%           categoryDecisionArr = sum(categoryArrs(nodes(1,1),:),1);
            % Mark previously considered best categories. We're selecting from
            % them.
            if ~isempty(bestCategories)
                categoryDecisionArr(setdiff(categoryList, bestCategories)) = 0;
            end
            [~, newBestCategories] = find(categoryDecisionArr == max(categoryDecisionArr));
            if ~isempty(bestCategories)
                bestCategories = intersect(bestCategories, newBestCategories);
            else
                bestCategories = newBestCategories;
            end
            if numel(bestCategories) == 1
                categoryLabel = bestCategories;
                break;
            end 
        end
%       
    end
end