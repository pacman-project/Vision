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
function [vocabLevel] = learnChildPositions(vocabLevel, graphLevel, modes)
    numberOfNodes = numel(vocabLevel);
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
        edges = vocabLevel(vocabItr).adjInfo;
        meanChildrenPos = zeros(numel(children),2, 'single');
        meanChildrenCov= zeros(numel(children),4, 'single');
        
        % Calculate mean positions for children and save them.
        for childItr = 2:numel(children)
             relevantMode = modes(modes(:,1) == children(1) & modes(:,2) == children(childItr) & modes(:,3) == edges(childItr-1, 3), :);

            % Now, we have the data. We can measure the mean and std of a
            % gaussian.
            meanChildrenPos(childItr,:) = relevantMode(10:11);
            meanChildrenCov(childItr,:) = relevantMode(6:9);
        end
        meanLocation = mean(meanChildrenPos,1);
        vocabLevel(vocabItr).childrenPosMean = meanChildrenPos - repmat(meanLocation, size(meanChildrenPos,1), 1);
        vocabLevel(vocabItr).childrenPosCov = meanChildrenCov;
    end
end