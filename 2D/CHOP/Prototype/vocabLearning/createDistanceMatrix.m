%> Name: createDistanceMatrix
%>
%> Description: Create initial distance matrix for CHOP. 
%> This function defines the similarity of level 1 parts, where 0 is the
%> closest (itself) and 1 is farthest.
%>
%> @param options Program options.
%>
%> @retval simMat Distance matrix of level 1 parts.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 01.09.2014
function [ distMat ] = createDistanceMatrix( filters )
    % Initialize internal data structures.
    numberOfFilters = numel(filters);
    distMat = zeros(numel(filters));

    % Find similarity between each pair of filters.
    for filtItr = 1:(numberOfFilters-1)
        filter1 = filters{filtItr};
        for filtItr2 = (filtItr+1):numberOfFilters
            filter2 = filters{filtItr2};
            distance = sqrt(sum(sum(sum((filter1-filter2).^2))));
            distMat(filtItr, filtItr2) = distance;
            distMat(filtItr2, filtItr) = distance;
        end
    end
    distMat = distMat/max(max(distMat));
end