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
function [ distMat, nodeLogMin ] = createDistanceMatrix( filters, filterType, distType, deadFeatures, desiredLogRange )
    distMat = zeros(numel(filters), 'single');
    nodeLogMin = 0;
    if strcmp(distType, 'prob') && strcmp(filterType, 'gabor')
        % We perform a search for optimal sigma here.
        lowSigma = 0.1;
        highSigma = 0.2;
        
        while true
             sigma = (lowSigma + highSigma)/2;
             % Calculate a 1D gaussian distribution here.
             numberOfFilters = numel(filters);
             anglePerFeature = 1 / double(numberOfFilters);
             halfAnglePerFeature = anglePerFeature / 2;
             startAngles = (0:(numberOfFilters-1)) * anglePerFeature - halfAnglePerFeature;
             endAngles = anglePerFeature + startAngles;
             refPoint = floor(numberOfFilters/2) + 1;
             refAngle = 0.5;

             % Measure probabilities and log probabilities
             startProbs = normcdf(startAngles, refAngle, sigma);
             endProbs = normcdf(endAngles, refAngle, sigma);
             probs = endProbs - startProbs;
             log2Probs = abs(log(probs));
             log2Probs = circshift(log2Probs, [0, refPoint-1]);
             
             % Learn log range.
             currentLogRange = max(log2Probs) - min(log2Probs);
             if currentLogRange > desiredLogRange + 0.001
                  lowSigma = sigma;
             elseif currentLogRange < desiredLogRange - 0.001
                  highSigma = sigma;
             else
                  break;
             end
        end
        
        % Scale log2Probs in a reversible way. We're not losing any info
        % here! We know the range and the minimum value. Chill!
        nodeLogMin = min(log2Probs);
        log2Probs = (log2Probs - min(log2Probs)) / currentLogRange;
        
        % Finally, assing the log probabilities to the distance matrix.
        for filtItr = 1:numberOfFilters
            assignedLogProbs = circshift(log2Probs, [0, (filtItr-1)]);
            distMat(filtItr,:) = assignedLogProbs;
        end
       return; 
    end
    
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
        trueFilter = filter1;
        filter1 = abs(filter1);
        for dimItr = 1:size(filter1,3)
            zeroIdx = filter1(:,:,dimItr) < max(max(max(filter1))) * 0.1;
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
        newFilter = circshift(trueFilter, cogDiff);
        cogFilters(filtItr) = {newFilter};
    end
    
    % Find distance between each pair of filters (cog-normalized).
    for filtItr = 1:(numberOfFilters-1)
        filter1 = cogFilters{filtItr};
        filter1 = filter1/norm(filter1(:));
        for filtItr2 = (filtItr+1):numberOfFilters
            filter2 = cogFilters{filtItr2};
            filter2 = filter2/norm(filter2(:));
            distance = findDistance(filter1, filter2, distType);
            distMat(filtItr, filtItr2) = distance;
            distMat(filtItr2, filtItr) = distance;
        end
    end
    % Normalize distMat.
    distMat(distMat == -1) = max(max(distMat));
    newDistMat = distMat/max(max(distMat));
    
   % Probability for auto mode has not been implemented yet!
    warning(12, 'Probability type distance for auto filters have not been implemented yet! Switching to euclidean measure..');
    
    % If rank type distance is used, each node's distances to others is
    % sorted, and the ranks are entered as the new distance functions.
    if strcmp(distType, 'rank')
        % We decide the relative order of similarity by taking each
        % filter and ranking the rest of the filters based on similarity. 
        newDistMat = zeros(size(distMat));
        sortAssgnArr = 1:numberOfFilters;
        for filtItr = 1:numberOfFilters
            distances = distMat(filtItr,:);
            [~, rankings] = sort(distances, 'ascend');
            [~,assgnArr] = ismember(sortAssgnArr, rankings);
            newDistMat(filtItr, :) = newDistMat(filtItr, :) + assgnArr;
            newDistMat(:, filtItr) = newDistMat(:, filtItr) + assgnArr';
        end
        newDistMat = newDistMat - 2;
        
%         % If the distance type is probabilities, we measure probabilities
%         % of filters.
%         if strcmp(distType, 'prob')
%              % We perform a search for optimal sigma here.
%              lowSigma = 0.05;
%              highSigma = 0.3;
%         
%              %% Find an optimal sigma so that the log range here matches the log range of edge case.
%              while true
%                   sigma = (lowSigma + highSigma)/2;
%                   numberOfFilters = numel(filters);
%                   anglePerFeature = 1 / double(numberOfFilters);
%                   halfAnglePerFeature = anglePerFeature / 2;
%                   startAngles = (0:(numberOfFilters-1)) * anglePerFeature - halfAnglePerFeature;
%                   endAngles = anglePerFeature + startAngles;
%                   refPoint = 1;
%                   refAngle = 0;
% 
%                   % Measure probabilities and log probabilities
%                   startProbs = normcdf(startAngles, refAngle, sigma);
%                   endProbs = normcdf(endAngles, refAngle, sigma);
%                   probs = endProbs - startProbs;
%                   probs = probs /  sum(probs);
%                   log2Probs = abs(log(probs));
%                   log2Probs = circshift(log2Probs, [0, refPoint-1]);
%                   
%                   % Learn log range.
%                   currentLogRange = max(log2Probs) - min(log2Probs);
%                   if currentLogRange > desiredLogRange + 0.001
%                        lowSigma = sigma;
%                   elseif currentLogRange < desiredLogRange - 0.001
%                        highSigma = sigma;
%                   else
%                        break;
%                   end
%              end
%              
%              % Scale log2Probs in a reversible way. We're not losing any info
%              % here! We know the range and the minimum value. Chill!
%              nodeLogMin = min(log2Probs);
%              log2Probs = log2Probs - min(log2Probs) / currentLogRange;
% 
%              % Finally, assing the log probabilities to the distance matrix.
%              for filtItr = 1:numberOfFilters
%    %              assignedLogProbs = circshift(log2Probs, [0, (filtItr-1)]);
%    %              newDistMat(filtItr,:) = assignedLogProbs;
%              end
%         end
    end
    validFeatures = setdiff(1:size(newDistMat), deadFeatures);
    distMat = newDistMat/max(max(newDistMat(validFeatures, validFeatures)));
    distMat(distMat > 1) = 1;
    distMat = single(distMat);
end

function distance = findDistance(filter1, filter2, distType)
    if strcmp(distType, 'euc') || strcmp(distType, 'rank')
        distance = sqrt(sum(sum(sum((filter1-filter2).^2))));
    else
        % This is where the probabilities for the node base case is
        % implemented.
        square_d = 0;
        XtY = multiprod(multitransp(filter1), filter2);
        for i = 1 : size(filter1,3)
             cos_princ_angle = svd(XtY(:, :, i));
             % Two next instructions not necessary: the imaginary parts that
             % would appear if the cosines are not between -1 and 1 when
             % passed to the acos function would be very small, and would
             % thus vanish when the norm is taken.
             % cos_princ_angle = min(cos_princ_angle,  1);
             % cos_princ_angle = max(cos_princ_angle, -1);
             square_d = square_d + norm(acos(cos_princ_angle))^2;
        end
         distance = sqrt(square_d);
    end
end