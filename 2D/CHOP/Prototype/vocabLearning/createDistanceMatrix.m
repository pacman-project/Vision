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
function [ distMat, abstractDistMat, nodeProbArr ] = createDistanceMatrix( filters, filterImages, filterType, distType )
    distMat = zeros(numel(filters), 'single');
    minFilterValue = 0.05; % If pixel a < max(max(max(filter1))) * minFilterValue,
                           % a is assigned 0 in distance calculations.
    gaborOrThr = 0;
    autoOrThr = 0.2;
    sigma = 0.005;
    cogFilters = cell(numel(filters),1);
    numberOfFilters = numel(filters);
    filterSize = [size(filters{1},1), size(filters{1},2)];
    binaryMask = true(filterSize);
    cog = zeros(1,2);
    trueCenter = 1 + (filterSize-1)/2;
    nodeProbArr = cell(numberOfFilters,1);
    
    % Find center of gravity for each filter, and center it around its
    % cog to estimate correct distance between pairs of filters.
    for filtItr = 1:numel(filters)
        filter1 = filters{filtItr};
        filterImage = filterImages{filtItr};
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
    %% Finally, we assign node substitution probabilities.
    if strcmp(filterType, 'gabor') 
         % Learn circular probabilities to assign for gabor choices.
         dataPoints = 1:(numberOfFilters+1);
         meanPoint = ceil((numberOfFilters+1)/2) / (numberOfFilters+2);
         startPoints = ((dataPoints - 1/2) / (numberOfFilters+2))';
         endPoints = ((dataPoints + 1/2) / (numberOfFilters+2))';
         startProbs = mvncdf(startPoints, meanPoint, sigma);
         endProbs = mvncdf(endPoints, meanPoint, sigma);
         pointProbs = endProbs - startProbs;
         
         % Finally, form a node probability array.
         for filtItr = 1:numberOfFilters
              entries = find(newDistMat(filtItr,:) <= gaborOrThr);
              realEntries = entries';
              entries = entries + round(numberOfFilters/2) - (filtItr - 1);
              entries = rem(entries + numberOfFilters, numberOfFilters);
              entries(entries==0) = numberOfFilters;
              probs = pointProbs(entries);
              
              % Normalize probabilities and save them.
              probs = probs/sum(probs);
              nodeProbArr{filtItr} = [realEntries, probs];
         end
    else
        for filtItr = 1:numberOfFilters
              % Rank closest filters, and assign probabilities based on
              % these ranks.
              closeFilters = find(newDistMat(filtItr,:) <= autoOrThr);
              distances = newDistMat(filtItr, closeFilters);
              [~, idx]=sort(distances,'ascend');
              sigma = 1/(numel(idx)*5);
              dataPoints = ((0:(numel(idx)-1))/numel(idx))';
              meanPoint = 0;
              startPoints = dataPoints - 1 / (2 * numel(idx));
              endPoints = startPoints + 1/numel(idx);
              startProbs = mvncdf(startPoints, meanPoint, sigma);
              endProbs = mvncdf(endPoints, meanPoint, sigma);
              pointProbs = endProbs - startProbs;
              
              % Normalize probabilities and save them.
              pointProbs = pointProbs/sum(pointProbs);
              nodeProbArr{filtItr} = [closeFilters(idx)', sort(pointProbs, 'descend')];
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