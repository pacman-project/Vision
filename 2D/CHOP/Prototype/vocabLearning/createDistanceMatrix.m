%> Name: createDistanceMatrix
%>
%> Description: Create initial distance matrix for CHOP. 
%> This function defines the similarity of level 1 parts, where 0 is the
%> closest (itself) and 1 is farthest.
%>
%> @param options Program options.
%>
%> @retval distMat Distance matrix of level 1 parts.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 01.09.2014
function [ distMat ] = createDistanceMatrix( filters )
    cogFilters = cell(numel(filters),1);
    numberOfFilters = numel(filters);
    distMat = zeros(numel(filters));
    filterSize = [size(filters{1},1), size(filters{1},2)];
    binaryMask = true(filterSize);
    cog = zeros(1,2);
    trueCenter = 1 + (filterSize-1)/2;
    % Find center of gravity for each filter, and center it around its
    % cog to estimate correct distance between pairs of filters.
    for filtItr = 1:numel(filters)
        filter1 = filters{filtItr};
        filter1(filter1<0) = filter1(filter1<0) * -1;
        filter1Gray = mean(filter1/ max(max(max(filter1))),3);
        % Assign small things to zero.
        minVal = max(max(filter1Gray)) * 0.2;
        filter1Gray(filter1Gray < max(max(filter1Gray)) * 0.2) = 0;
        measurements = regionprops(binaryMask, filter1Gray, 'WeightedCentroid');
        cog([2 1]) = measurements(1).WeightedCentroid;
        cogDiff = round(trueCenter - cog);
        newFilter = circshift(filters{filtItr}, cogDiff);
        newFilter(newFilter < minVal & newFilter > -minVal) = 0;
        cogFilters(filtItr) = {newFilter};
    end
    % Find distance between each pair of filters (cog-normalized).
    for filtItr = 1:(numberOfFilters-1)
        filter1 = cogFilters{filtItr};
        for filtItr2 = (filtItr+1):numberOfFilters
            filter2 = cogFilters{filtItr2};
            distance = sqrt(sum(sum(sum((filter1-filter2).^2))));
            distMat(filtItr, filtItr2) = distance;
            distMat(filtItr2, filtItr) = distance;
        end
    end
    % Normalize distMat.
    distMat = distMat/max(max(distMat));
end