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
%> Ver 1.0 on 05.05.2014
%> Ver 1.1 on 12.08.2014 Adding "auto" type similarity matrix.
function [ simMat ] = createSimilarityMatrix( options )
    if strcmp(options.filterType, 'gabor') || strcmp(options.filterType, 'lhop')
        %% TODO: Change (CHANGE) the following similarity matrix definition with something more 'mathematical'. 
        % Oh god.
        simMat = [0 1/3 2/3 1 2/3 1/3; ...
                1/3 0 1/3 2/3 1 2/3; ...
                2/3 1/3 0 1/3 2/3 1; ... 
                1 2/3 1/3 0 1/3 2/3; ...
                2/3 1 2/3 1/3 0 1/3; ...
                1/3 2/3 1 2/3 1/3 0];
        % Oh god.. Shakespeare would be proud.
    elseif strcmp(options.filterType, 'auto')
        filters = options.filters;
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
    else
        simMat = [];
    end
end