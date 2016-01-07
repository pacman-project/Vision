%> Name: learnChildDistributions
%>
%> Description: This function learns label/position distributions of the sub-parts
%> of every part in vocabLevel, using the data in the object graphs.
%> 
%> @param vocabLevel Current vocabulary level.
%> @param graphLevel Current graph level.
%> @param previousLevel Previous graph level.
%> 
%> @retval vocabLevel The vocabulary level updated with sub-part
%> distributions.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 07.07.2015
function [vocabLevel] = learnChildDistributions(vocabLevel, graphLevel, previousLevel)
    numberOfNodes = numel(vocabLevel);
    labelIds = [graphLevel.labelId];
    curPrecisePos = cat(1, graphLevel.precisePosition);
    prevRealLabelIds = [previousLevel.realLabelId];
    prevPrecisePos = cat(1, previousLevel.precisePosition);
    
    % Second, we go through each node, and collect statistics.
    for vocabItr = 1:numberOfNodes
        children = vocabLevel(vocabItr).children;
        meanChildrenPos = zeros(numel(children),2, 'single');
        
        % Get data. 
        % TODO: Consider mapping here!
        relevantIdx = labelIds == vocabItr;
        instances = graphLevel(relevantIdx);
        instanceChildren = cat(1, instances.children);
        for instanceItr = 1:size(instanceChildren,2)
             relevantChildren = instanceChildren(:, instanceItr);
             childLabel = mode(double(prevRealLabelIds(relevantChildren)));
             children(instanceItr) = childLabel;
             meanChildrenPos(instanceItr,:) = mean(prevPrecisePos(relevantChildren, :) - curPrecisePos(relevantIdx, :),1);
        end
        vocabLevel(vocabItr).childrenPosMean = meanChildrenPos;
        vocabLevel(vocabItr).realChildren = children;
    end
end