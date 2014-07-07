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
    maxLevels = max(exportArr(:,4));
    categoryLabel = -1;
    for levelItr = maxLevels:-1:1
        vocabLevel = vocabulary{levelItr};
        categoryArrs = cat(1, vocabLevel.categoryArr);
        if isempty(categoryArrs) 
           break; 
        end
        nodes = exportArr(exportArr(:,4) == levelItr,:);
        categoryDecisionArr = sum(categoryArrs(nodes(:,1),:),1) / size(nodes,1);
        [~, bestCategories] = find(categoryDecisionArr == max(categoryDecisionArr));
        
        if numel(bestCategories) == 1
            categoryLabel = bestCategories;
            break;
        end
    end
end