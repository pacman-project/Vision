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
%>  
%> @retval graphLevel 
%>
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 04.02.2014
function graphLevel = fillBasicInfo(previousLevel, graphLevel, leafNodes)
    for newNodeItr = 1:numel(graphLevel)
        graphLevel(newNodeItr).imageId = previousLevel(graphLevel(newNodeItr).children(1)).imageId;
        numberOfChildren = numel(graphLevel(newNodeItr).children);
        nodeLeafNodes = cell(numberOfChildren,1);
        for childItr = 1:numel(graphLevel(newNodeItr).children)
           nodeLeafNodes(childItr) = {previousLevel(graphLevel(newNodeItr).children(childItr)).leafNodes}; 
        end
        nodeLeafNodes = unique([nodeLeafNodes{:}]);
        graphLevel(newNodeItr).position = round(sum(cat(1, leafNodes{nodeLeafNodes,2}),1) / numel(nodeLeafNodes));
        graphLevel(newNodeItr).leafNodes = nodeLeafNodes;
    end

    % Rearrange graph level so it is sorted by image id.
    arrayToSort = [[graphLevel.imageId]', [graphLevel.labelId]'];
    [~, sortedIdx] = sortrows(arrayToSort);
    graphLevel = graphLevel(sortedIdx);
end