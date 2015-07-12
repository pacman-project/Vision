%> Name: learnChildPositions
%>
%> Description: This function learns the mean and standard deviation of the
%> gaussian distribution imposed by each child of vocabulary nodes. For each
%> child, relative positioning data is collected from the realizations of
%> this node. Then, a gaussian is fit to the relation of every child. 
%>
%> @param vocabLevel Current vocabulary level.
%> @param graphLevel Current graph level.
%> @param previousLevel Previous graph level.
%> 
%> @retval vocabLevel The vocabulary level updated with children mean
%> positions/standard deviations. 
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 07.07.2015
function [vocabLevel] = learnChildPositions(vocabLevel, graphLevel, previousLevel)
    numberOfNodes = numel(vocabLevel);
    previousLevelPos = single(cat(1, previousLevel.position));
    
    %% For every vocabulary node, go through all the graph nodes and learn statistics of children.
    % First, we obtain instances for every abstract node.
    graphNodeArr = cell(numberOfNodes,1);
    labelIds = [graphLevel.labelId];
    for vocabItr = 1:numberOfNodes
        graphNodeArr(vocabItr) = {graphLevel(labelIds == vocabItr)};
    end
    
    % Second, we go through each node, and collect statistics.
    for vocabItr = 1:numberOfNodes
        children = vocabLevel(vocabItr).children;
        instances = graphNodeArr{vocabItr};
        instanceChildren = cat(1, instances.children);
        instancePrecisePositions = cat(1, instances.precisePosition);
        instanceMappings = cat(1, instances.mapping);
        meanChildrenPos = zeros(numel(children),2, 'single');
        meanChildrenStd = zeros(numel(children),1, 'single');
        
        % Calculate mean positions for children and save them.
        for childItr = 1:numel(children)
            assignedChildren = instanceChildren(instanceMappings==childItr);
            assignedChildrenPos = instancePrecisePositions - previousLevelPos(assignedChildren, :);
            
            % Now, we have the data. We can measure the mean and std of a
            % gaussian.
            meanChildrenPos(childItr,:) = mean(assignedChildrenPos, 1);
            
            % Get std. 
            relativeChildrenPos = assignedChildrenPos - ...
                repmat(meanChildrenPos(childItr,:), size(assignedChildrenPos,1), 1);
            meanChildrenStd(childItr,:) = sqrt(sum(sum((relativeChildrenPos.^2), 2)) ...
                / size(assignedChildrenPos,1));
        end
        vocabLevel(vocabItr).childrenPosMean = meanChildrenPos;
        vocabLevel(vocabItr).childrenPosStd = meanChildrenStd;
    end
end