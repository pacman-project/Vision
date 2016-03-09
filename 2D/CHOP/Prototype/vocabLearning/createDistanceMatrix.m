%> Name: createDistanceMatrix
%>
%> Description: Create initial distance matrix for CHOP. 
%> This function defines the similarity of level 1 parts, where 0 is the
%> closest (itself) and 1 is farthest.
%>
%> @param filters The masks of level 1 filters of the hierarchy in a cell %
%> array.
%> @param filterType The type of the filters could be either gabor or auto.
%> @param distType If 'euc', Euclidean distances of two masks (after their
%> centers of gravity are aligned), normalized by the max distance is used as
%> pairwise node substitution cost. If 'rank', for each node, other nodes
%> are ordered by their pairwise distances to this node, and their ranks are
%> used as the distance function. If probability, a von mises distribution
%> is used to calculate node replacement probabilities.
%> @param deadFeatures If non-empty, we eliminate the filters who have ids 
%> in this array from our calculations. 
%>
%> @retval distMat Distance matrix of level 1 parts.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 01.09.2014
%> Comments updated on 23.06.2015
function [ distMat, abstractDistMat ] = createDistanceMatrix( filters, filterType)
    distMat = zeros(numel(filters), 'single');
    minFilterValue = 0.05; % If pixel a < max(max(max(filter1))) * minFilterValue,
                           % a is assigned 0 in distance calculations.
    gaborOrThr = 0;
    autoOrThr = 0.25;
    cogFilters = cell(numel(filters),1);
    numberOfFilters = numel(filters);
    filterSize = [size(filters{1},1), size(filters{1},2)];
    binaryMask = true(filterSize);
    cog = zeros(1,2);
    trueCenter = 1 + (filterSize-1)/2;
    
    % Find center of gravity for each filter, and center it around its
    % cog to estimate correct distance between pairs of filters.
    for filtItr = 1:numel(filters)
        filter1 = filters{filtItr};
 %       filterImage = filterImages{filtItr};
        filterImage = filter1;
        filterImage = (filterImage - min(min(min(filterImage)))) / (max(max(max(filterImage))) - min(min(min(filterImage))));
        trueFilter = filter1;
        filter1 = abs(filter1);
        for dimItr = 1:size(filter1,3)
            zeroIdx = filter1(:,:,dimItr) < max(max(max(filter1))) * minFilterValue;
            channelImg = filter1(:,:,dimItr);
            channelImg(zeroIdx) = 0;
            filter1(:,:,dimItr) = channelImg;
            channelImg = trueFilter(:,:,dimItr);
            channelImg(zeroIdx) = 0;
            trueFilter(:,:,dimItr) = channelImg;
        end
        filter1Gray = mean(filter1/ max(max(max(filter1))),3);
        measurements = regionprops(binaryMask, filter1Gray, 'WeightedCentroid');
        cog([2 1]) = measurements(1).WeightedCentroid;
        cogDiff = round(trueCenter - cog);
        newFilter = circshift(double(filterImage), cogDiff);
        cogFilters(filtItr) = {newFilter};
    end
    
    % Find distance between each pair of filters (cog-normalized).
    for filtItr = 1:(numberOfFilters-1)
        filter1 = cogFilters{filtItr};
        for filtItr2 = (filtItr+1):numberOfFilters
            filter2 = cogFilters{filtItr2};
            distance = findDistance(filter1, filter2);
            distMat(filtItr, filtItr2) = distance;
            distMat(filtItr2, filtItr) = distance;
        end
    end
    % Normalize distMat.
    distMat(distMat == -1) = max(max(distMat));
    newDistMat = distMat/max(max(distMat));
    
    % Calculate OR nodes for first level filters. We are pretty flexible
    % in the first level.
    if strcmp(filterType, 'gabor') 
        % We decide the relative order of similarity by taking each
        % filter and ranking the rest of the filters based on similarity. 
        newDistMat = zeros(size(distMat));
        firstRow = [0:floor(numberOfFilters/2), floor((numberOfFilters-1)/2):-1:1];
        firstRow = firstRow/max(firstRow);
        for filtItr = 1:numberOfFilters
             newDistMat(filtItr,:) = circshift(firstRow', (filtItr-1))'; 
        end
    end
    
    %% Save valid features and normalize the distance matrix.
%    validFeatures = setdiff(1:size(newDistMat), deadFeatures);
%    distMat = newDistMat(validFeatures, validFeatures);
    distMat = single(newDistMat);
    abstractDistMat = distMat;
    if strcmp(filterType, 'gabor') 
         abstractDistMat(abstractDistMat > gaborOrThr) = inf;
    else
         abstractDistMat(abstractDistMat > autoOrThr) = inf;
    end
    abstractDistMat(abstractDistMat < inf) = 0;
end

function distance = findDistance(filter1, filter2)
    distance = sqrt(sum(sum(sum((filter1-filter2).^2))));
end