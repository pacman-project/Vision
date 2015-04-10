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
function [ categoryLabel, decisionLevel ] = getCategoryLabel(vocabulary, exportArr, activationArr, minLevels, maxLevels)
%    maxLevels = max(exportArr(:,4));
    categoryLabel = -1;
    numberOfCategories = numel(vocabulary{1}(1).categoryArr);
    if ~isempty(activationArr)
        activationArr = repmat(activationArr, 1, numberOfCategories);
    end
    maxLevels = min(numel(vocabulary), maxLevels);
    for levelItr = maxLevels:-1:minLevels
        categoryDecisionArr = zeros(1, numberOfCategories, 'single');
        vocabLevel = vocabulary{levelItr};
        categoryArrs = double(cat(1, vocabLevel.categoryArr)).^2;
        if isempty(categoryArrs) 
           break; 
        end
        nodes = exportArr(exportArr(:,4) == levelItr,:);
        if ~isempty(activationArr)
             nodeActivations = activationArr(exportArr(:,4) == levelItr,:);
%             % Here, we set the maximum activation as 1, and others as zero.
%             if size(nodeActivations,1) > 1
%                 [maxVal, ~] = max(max(nodeActivations));
%                 nodeActivations(nodeActivations<(maxVal-0.000001)) = 0;
%             end
        else
            nodeActivations = ones(size(exportArr,1),1);
        end
        if size(nodes,1)>0
            categoryDecisionArr = categoryDecisionArr + sum(categoryArrs(nodes(:,1),:) .* nodeActivations, 1);
            [~, newBestCategories] = find(categoryDecisionArr == max(categoryDecisionArr));
            if numel(newBestCategories) == 1
                categoryLabel = newBestCategories;
                break;
            end 
        end
    end
    decisionLevel = levelItr;
end