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
function graphLevel = fillBasicInfo(previousLevel, graphLevel, leafNodes, numberOfThreads)
    numberOfNodes = numel(graphLevel);
    nodeSets = repmat(ceil(numberOfNodes/numberOfThreads), numberOfThreads,1);
    remNodes = rem(numberOfNodes, numberOfThreads);
    if remNodes > 0
        nodeSets(end) = nodeSets(end) + remNodes - numberOfThreads;
    end
    if numel(graphLevel) >= ceil(numberOfNodes/numberOfThreads) * numberOfThreads
        nodeSets = mat2cell(graphLevel, 1, nodeSets);
    else
        nodeSets = {graphLevel};
    end
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
            subLevel(newNodeItr).position = int32(round(sum(cat(1, leafNodes(nodeLeafNodes,2:3))) / numel(nodeLeafNodes)));
            subLevel(newNodeItr).leafNodes = nodeLeafNodes;
        end
        nodeSets(setItr) = {subLevel};
    end
    graphLevel = cat(2, nodeSets{:});

    % Rearrange graph level so it is sorted by image id.
    arrayToSort = [[graphLevel.imageId]', [graphLevel.labelId]'];
    [~, sortedIdx] = sortrows(arrayToSort);
    graphLevel = graphLevel(sortedIdx);
end