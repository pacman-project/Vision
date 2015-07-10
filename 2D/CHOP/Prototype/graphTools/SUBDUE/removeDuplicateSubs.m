%> Name: removeDuplicateSubs
%>
%> Description: Given the child subs in childSubArr, this function removes
%> duplicate substructures from childSubArr, and combines the instances of
%> matching subs. If an instance can be parsed in different ways (i.e.
%> matches multiple duplicate subs), the minimum matching cost of matching is
%> saved. 
%> 
%> @param childSubArr A list of children substructures.
%>
%> @retval childSubArr The list of unique substructures.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 01.03.2015 (Converted to a standalone function, added instance augmentation)
function [childSubArr] = removeDuplicateSubs(childSubArr)
    if size(childSubArr(1).edges,1) > 1
        % Eliminate duplicate subs in final array.
        childSubEdges = cell(1, numel(childSubArr));
        for subItr = 1:numel(childSubArr)
            childSubEdges(subItr) = {sortrows(childSubArr(subItr).edges)};
        end
        subCenters = {childSubArr.centerId};
        vocabDescriptors = cellfun(@(x,y) num2str([x; y(:)]'), subCenters, childSubEdges, 'UniformOutput', false);
        [~, validChildrenIdx, ~] = unique(vocabDescriptors, 'stable');
        childSubArr = childSubArr(validChildrenIdx);
    end
end