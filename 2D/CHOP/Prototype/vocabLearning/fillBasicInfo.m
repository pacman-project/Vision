%> Name: fillBasicInfo
%>
%> Description: Given current level and previous level (both object graphs), 
%> fill in the current level's necessary basic info, including image id,
%> position, and leafNodes. In the end, the output is sorted based on imageId
%> and labelId.
%>
%> @param previousLevel 
%> @param graphLevel 
%> @param leafNodes 
%> @param numberOfThreads
%>  
%> @retval graphLevel 
%>
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 04.02.2014
function graphLevel = fillBasicInfo(previousLevel, graphLevel, ~, numberOfThreads)
    numberOfNodes = numel(graphLevel);
    nodeSets = repmat(ceil(numberOfNodes/numberOfThreads), numberOfThreads,1);
    setCountDiff = sum(nodeSets) - numberOfNodes;
    nodeSets((end-setCountDiff+1):end, :) = nodeSets((end-setCountDiff+1):end, :) - 1;
    nodeSets = nodeSets(nodeSets>0);
    nodeSets = mat2cell(graphLevel, 1, nodeSets);
    parfor setItr = 1:numel(nodeSets)
        subLevel = nodeSets{setItr};
        for newNodeItr = 1:numel(subLevel)
            nodeChildren = subLevel(newNodeItr).children;
            numberOfChildren = numel(nodeChildren);
            subLevel(newNodeItr).imageId = previousLevel(nodeChildren(1)).imageId;
            nodeLeafNodes = cell(numberOfChildren,1);
            for childItr = 1:numberOfChildren
               nodeLeafNodes(childItr) = {previousLevel(nodeChildren(childItr)).leafNodes}; 
            end
            nodeLeafNodes = unique([nodeLeafNodes{:}]);
            precisePosition = sum(cat(1, previousLevel(nodeChildren).precisePosition),1) ...
                / numberOfChildren;
            subLevel(newNodeItr).precisePosition = precisePosition;
            subLevel(newNodeItr).position = int32(round(mean(cat(1, previousLevel(nodeChildren).position),1)));
            subLevel(newNodeItr).leafNodes = nodeLeafNodes;
        end
        nodeSets(setItr) = {subLevel};
    end
    graphLevel = cat(2, nodeSets{:});

    % Rearrange graph level so it is sorted by image id.
    arrayToSort = [[graphLevel.imageId]', [graphLevel.labelId]', cat(1, graphLevel.position)];
    [~, sortedIdx] = sortrows(arrayToSort);
    graphLevel = graphLevel(sortedIdx);
end