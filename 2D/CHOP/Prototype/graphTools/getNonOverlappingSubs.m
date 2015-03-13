%> Name: getNonOverlappingSubs
%>
%> Description: Given the extended subs in childSubs, this function tries
%> to eliminate the compositions which match (within the limits specified by
%> adaptiveThreshold) to more frequent subs in childSubs. The purpose is to
%> lower the computational complexity, without giving away too much in
%> precision.
%> 
%> @param childSubs A list of substructures.
%> @param nodeDistanceMatrix The distance matrix of the nodes in the
%> previous layer. 
%> @param edgeDistanceMatrix The edge distance matrix which shows how
%> geometrically distant the edge types are.
%> @param adaptiveThreshold The threshold that specifies elasticity. It is
%> calculated using the size of the subs.
%> @param singlePrecision Precision of the numeric values.
%>
%> @retval childSubs The list of unique substructures.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 25.02.2015 (Converted to a standalone function)
function [childSubs, validSubs] = getNonOverlappingSubs(childSubs, nodeDistanceMatrix, edgeDistanceMatrix, adaptiveThreshold, singlePrecision)
    validSubs = ones(numel(childSubs),1) > 0;
    if adaptiveThreshold > singlePrecision
        numberOfChildSubs = numel(childSubs);
        
        % If they have not been evaluated yet, we order them by occurence
        % counts.
        if isempty(childSubs(1).mdlScore)
            freqScores = zeros(numberOfChildSubs,1);
            for subItr = 1:numberOfChildSubs
                freqScores(subItr) = size(childSubs(subItr).instanceCenterIdx,1);
            end
            [~, sortIdx] = sort(freqScores, 'descend');
            childSubs = childSubs(sortIdx);
        end
        
        % Find valid subs.
        edgeNodePairs = cat(1, childSubs.edges);
        numberOfEdgesPerSub = size(childSubs(1).edges,1);
        numberOfChildSubs = numel(childSubs);
        edgeNodePairs = edgeNodePairs(numberOfEdgesPerSub:numberOfEdgesPerSub:(numberOfChildSubs*numberOfEdgesPerSub), :);
        for childItr = 1:(numberOfChildSubs-1)
            if validSubs(childItr)
                edgeNodePair1 = edgeNodePairs(childItr,:);
                for childItr2 = (childItr+1):numberOfChildSubs
                    if validSubs(childItr2)
                        edgeNodePair2 = edgeNodePairs(childItr2,:);
                        if edgeDistanceMatrix(edgeNodePair1(1), edgeNodePair2(1)) + ...
                                nodeDistanceMatrix(edgeNodePair1(2), edgeNodePair2(2)) < adaptiveThreshold
                            validSubs(childItr2) = 0;
                        end
                    end
                end
            end
        end
        childSubs = childSubs(validSubs);
    end
end