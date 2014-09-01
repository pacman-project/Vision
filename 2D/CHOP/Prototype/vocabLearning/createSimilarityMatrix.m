%> Name: createSimilarityMatrix
%>
%> Description: Create initial distance matrix for CHOP. 
%> This function defines the similarity of level 1 parts, where 0 is the
%> closest (itself) and 1 is farthest.
%> TODO: Change name to createDistanceMatrix.
%>
%> @param options Program options.
%>
%> @retval simMat Distance matrix of level 1 parts.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 01.09.2014
function [ simMat ] = createSimilarityMatrix( filters )
    cogFilters = cell(numel(filters),1);
    numberOfFilters = numel(filters);
    simMat = zeros(numel(filters));
    filterSize = [size(filters{1},1), size(filters{1},2)];
    binaryMask = true(filterSize);
    cogs = zeros(numel(filters),2);
    trueCenter = round(filterSize/2);

    % Find center of gravity for each filter, and center it around its
    % cog to estimate correct distance between pairs of filters.
    for filtItr = 1:numel(filters)
        filter1 = filters{filtItr};

        if size(filter1,3) > 1
            filter1Gray = mean(filter1,3);
        else
            filter1Gray = filter1;
        end
        filter1Gray = (filter1Gray-min(min(filter1Gray)))/ (max(max(filter1Gray)) - min(min(filter1Gray)));
        measurements = regionprops(binaryMask, filter1Gray, 'WeightedCentroid');
        cogs(filtItr,[2 1]) = round(measurements(1).WeightedCentroid);
        halfDims = min((cogs(filtItr,:)-1), filterSize - cogs(filtItr,:));
        newFilter = zeros(size(filter1), class(filter1));
        newFilter((trueCenter(1)-halfDims(1)):(trueCenter(1)+halfDims(1)), ...
            (trueCenter(2)-halfDims(2)):(trueCenter(2)+halfDims(2)), :) = ...
            filter1((cogs(filtItr,1)-halfDims(1)):(cogs(filtItr,1)+halfDims(1)), ...
            (cogs(filtItr,2)-halfDims(2)):(cogs(filtItr,2)+halfDims(2)), :);
        cogFilters(filtItr) = {newFilter};
    end

    % Find similarity between each pair of filters (cog-normalized).
    for filtItr = 1:(numberOfFilters-1)
        filter1 = cogFilters{filtItr};
        for filtItr2 = (filtItr+1):numberOfFilters
            filter2 = cogFilters{filtItr2};
            distance = sqrt(sum(sum(sum((filter1-filter2).^2))));
            simMat(filtItr, filtItr2) = distance;
            simMat(filtItr2, filtItr) = distance;
        end
    end
    simMat = simMat/max(max(simMat));
end