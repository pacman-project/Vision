%> Name: getDisjointSubs
%>
%> Description: Given a set of parts in bestSubs, this function matches
%> parts with common ancestors. If two parts only differ in the last added 
%> edge/node pair, a pairwise distance between them is calculated. If their
%> pairwise distance is smaller than a threshold, they match. In this case,
%> the one with lower number of leaf nodes is eliminated. 
%> 
%> @param bestSubs A set of substructures evaluated on the training set.
%> @param subLeafNodes The leaf node list for each substructure.
%> @param threshold Matching threshold. 
%>
%> @retval validSubIdx Linear indices of valid (kept) subs.
%>
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 23.02.2015
function [ validSubIdx ] = getDisjointSubs( bestSubs, ...
    nodeDistanceMatrix, edgeDistanceMatrix, singlePrecision, threshold )

    numberOfBestSubs = numel(bestSubs);
    multipleNodeSubIdx = ones(numberOfBestSubs,1)>0;
    validSubIdx = zeros(numberOfBestSubs,1)>0;

    % Find subs with common ancestors.
    childSubEdges = cell(1, numberOfBestSubs);
    for subItr = 1:numberOfBestSubs
        if size(bestSubs(subItr).edges,1) > 1
            childSubEdges(subItr) = {bestSubs(subItr).edges(1:(end-1),:)};
        elseif isempty(bestSubs(subItr).edges)
            multipleNodeSubIdx(subItr) = 0;
        end
    end
    validSubIdx(~multipleNodeSubIdx) = 1;
    subCenters = {bestSubs.centerId};
    subCenters(~multipleNodeSubIdx) = {[]};
    vocabDescriptors = cellfun(@(x,y) num2str([x; y(:)]'), subCenters, childSubEdges, 'UniformOutput', false);
    [subRoots, validChildrenIdx, IC] = unique(vocabDescriptors, 'stable');
    
    % Count number of occurences for each sub root. For any one which has
    % been seen more than once, we will run an elimination of overlapping
    % subs.
    occurenceArr = hist(IC, 1:numel(validChildrenIdx));
    subRootsValid = occurenceArr == 1;
    validSubIdx(validChildrenIdx(subRootsValid)) = 1;
    subRootsToTest = find(occurenceArr > 1 & cellfun(@(x) ~isempty(x), subRoots));
    
    % Go through every possible root, and eliminate overlapping subs.
    for subRootItr = subRootsToTest
        subsToTest = IC == subRootItr;
        subListToTest = bestSubs(subsToTest);
        adaptiveThreshold = (2 * size(subListToTest(1).edges,1) + 1) * threshold + singlePrecision;
        [~, remainingSubIdx] = getNonOverlappingSubs(bestSubs(subsToTest), ...
            nodeDistanceMatrix, edgeDistanceMatrix, adaptiveThreshold, singlePrecision);
        subsToTestIdx = find(subsToTest);
        remainingSubIdx = subsToTestIdx(remainingSubIdx);
        validSubIdx(remainingSubIdx) = 1;
    end
end

