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
 %   overlapThr = 0.9;

    if size(childSubArr(1).edges,1) > 0
        % Eliminate duplicate subs in final array.
        childSubEdges = cell(1, numel(childSubArr));
        for subItr = 1:numel(childSubArr)
            childSubEdges(subItr) = {sortrows(childSubArr(subItr).edges)};
        end
        subCenters = {childSubArr.centerId};
        vocabDescriptors = cellfun(@(x,y) num2str([x; y(:)]'), subCenters, childSubEdges, 'UniformOutput', false);
        [~, validChildrenIdx, ~] = unique(vocabDescriptors, 'stable');
        childSubArr = childSubArr(validChildrenIdx);
        
%         % If finalListFlag is on, we perform another set of operations to
%         % ensure mostly overlapping subs do not make the final list.
%         if finalListFlag
%            display(['[SUBDUE] Removing highly overlapping subs from final list (Threshold is set to ' num2str(overlapThr) ').']);
%            subChildren = arrayfun(@(x) [x.centerId x.edges(:,2)'], childSubArr, 'UniformOutput', false);
%            subChildren = sort(cat(1, subChildren{:}), 2);
%            uniqueChildren = {childSubArr.instanceChildren};
%            uniqueChildren = cellfun(@(x) unique(x), uniqueChildren, 'UniformOutput', false);
%            uniqueChildrenCounts = cellfun(@(x) numel(x), uniqueChildren)';
%            validSubs = ones(size(subChildren,1),1) > 0;
%            for subItr = 1:numel(childSubArr) 
%                if validSubs(subItr)
%                    subsToMatch = validSubs;
%                    subsToMatch = subsToMatch & uniqueChildrenCounts >= uniqueChildrenCounts(subItr)*overlapThr & ...
%                        uniqueChildrenCounts <= uniqueChildrenCounts(subItr)*(1/overlapThr);
%                    subsToMatch(1:subItr) = 0;
%                    subsToMatchIds = find(subsToMatch);
%                    candidateMatchingSubs = ismember(subChildren(subsToMatch, :), subChildren(subItr,:), 'rows');
%                    if nnz(candidateMatchingSubs) > 0
%                        subsToCheck = subsToMatchIds(candidateMatchingSubs);
%                        checkSubArr = cellfun(@(x) numel(fastintersect(uniqueChildren{subItr}, x)) <= overlapThr * numel(uniqueChildren{subItr}), uniqueChildren(subsToCheck));
%                        validSubs(subsToCheck) = checkSubArr;
%                    end
%                end
%            end
%            childSubArr = childSubArr(validSubs);
%         end
    end
end