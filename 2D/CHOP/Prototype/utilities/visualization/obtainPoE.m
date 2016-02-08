%> Name: obtainPoE
%>
%> Description: Given level 1 nodes' projections, we calculate a Product of
%> Experts (PoE) agreement over image pixels.
%>
%> @param level1Nodes exported nodes in the format of 
%> [labelId, centerPosX, centerPosY].
%> @param imgSize Image size in pixels.
%> @param options Program options.
%> @param fastFlag Fast processing flag.
%> @param muImgInitial Initial predictions (modal) for each pixel, if there are any. 
%> Used to speed up calculation process.  
%> @param varImgInitial Variance of the predictions for each pixel, if any.
%>  
%> @retval muImg Final predictions.
%> @retval varImg Variance of final predictions.
%>
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 07.02.2016
function [ muImg, varImg ] = obtainPoE( level1Nodes, imgSize, options, fastFlag, muImgInitial, varImgInitial)
     imgSize = double(imgSize);
     filterSize = size(options.filters{1},1);
     halfFilterSize = round((size(options.filters{1},1) - 1) / 2);
     maxInfluenceRadius = ceil(halfFilterSize * sqrt(2));
     
     % If initial mu/var images are provided, we're only supposed to update
     % what we have.
     if nargin == 4
          muImg = zeros(imgSize);
          varImg = zeros(imgSize);
     else
          muImg = muImgInitial;
          varImg = varImgInitial;
     end
     
     % Eliminate unnecessary nodes.
     level1Nodes = level1Nodes(level1Nodes(:,2) > 1 & level1Nodes(:,2) < imgSize(1) + 1 & ...
           level1Nodes(:,3) > 1 & level1Nodes(:,3) < imgSize(2) + 1, :);

     % Create overlap array.
     overlapMatrix = zeros(imgSize) > 0;
     idx = sub2ind(imgSize, level1Nodes(:,2), level1Nodes(:,3));
     overlapMatrix(idx) = 1;
     maxOverlapMatrix = imdilate(overlapMatrix, strel('disk', maxInfluenceRadius));
     overlapMatrix = imdilate(overlapMatrix, strel('square', filterSize));
     
     % Find sigma values for every distance.
     maxDist = ceil(sqrt(imgSize(1)^2 + imgSize(2)^2)) + 1;
     distVals = 0:1/(maxInfluenceRadius):1;
     gaussVals = normpdf(distVals, 0, 0.4);
     allSigmaVals = 1./gaussVals;
     allSigmaVals = allSigmaVals(1:(maxInfluenceRadius+1))';
     
     % Calculate remaining program variables.
     posArr = double(level1Nodes(:,2:3));
     gaborIdArr = double(level1Nodes(:,1));
     numberOfGabors = size(level1Nodes,1);
     if size(level1Nodes,2) > 3
          activationArr = level1Nodes(:,4);
     else
          activationArr = ones(numberOfGabors,1);
     end
     
     filters = options.filters;
     filters = cellfun(@(x) (x - min(min(x))) / (max(max(x)) - min(min(x))), filters, 'UniformOutput', false);
     firstHalf = ceil(options.gaborFilterSize/2) - 1;
     secHalf = options.gaborFilterSize - (firstHalf+1);
     
     % Calculate pixel-level likelihoods.
     for itr1 = 1:imgSize(1)
          for itr2 = 1:imgSize(2)
               if (fastFlag && ~overlapMatrix(itr1, itr2)) || (~fastFlag && ~maxOverlapMatrix(itr1, itr2))
                    continue;
               end
               
               location = [itr1, itr2];
               distances = repmat(location, numberOfGabors, 1) - posArr;
               actualDistances = floor(sqrt(sum(distances.^2,2))) + 1;
               actualDistances(actualDistances > maxDist) = maxDist;
               invalidIdx = actualDistances > maxInfluenceRadius;
               sigmaVals = zeros(numel(actualDistances),1);
               sigmaVals(~invalidIdx) = allSigmaVals(actualDistances(~invalidIdx));
               overlappingIds = distances(:,1) >(-firstHalf-1) & distances(:,1) < (secHalf+1) & distances(:,2) >(-firstHalf-1) & distances(:,2) < (secHalf+1);

               % Keep track of overlapping gabors filters, and use them in
               % calculations.
               overlappingIdx = find(~invalidIdx);
               
               % No experts? Move on. We're confident.
               if isempty(overlappingIdx)
                    aggMu = 0;
                    aggSigma = 0;
               else 
                    %% Get experts.
                    % Allocate space for predictions.
                    predictions = zeros(numel(overlappingIdx),2);
                    for expItr = 1:numel(overlappingIdx)
                         expId = overlappingIdx(expItr);
                         % Get the mu of the new distribution.
                         if overlappingIds(expId) == 1
                              filterId = gaborIdArr(expId);
                              filterVals = filters{filterId};
                              validPos = distances(expId,:) + [firstHalf, secHalf] + 1;
                              predictions(expItr,1) = filterVals(validPos(1), validPos(2)) * activationArr(expId);   
                         else
                              predictions(expItr,1) = muImg(itr1, itr2);
                         end
                         predictions(expItr,2) = sigmaVals(expId);
                    end
                    
                    % Obtain earlier predictions.
                    if varImg(itr1, itr2) > 0
                         predictions = [muImg(itr1, itr2), varImg(itr1, itr2); predictions]; %#ok<AGROW>
                    end

                    %% Calculate product of experts.
                    aggMu = predictions(1,1);
                    aggSigma = predictions(1,2);
                    for expItr = 2:size(predictions,1)
                         newMu = predictions(expItr,1);
                         newSigma = predictions(expItr,2);
                         % Take product of gaussians.
                         aggMu = (aggMu * (newSigma^2) + newMu * (aggSigma^2)) / (aggSigma^2 + newSigma^2);
                         aggSigma = sqrt(((aggSigma^2)*(newSigma^2))/(aggSigma^2+newSigma^2));
                    end
               end
               
               % Assign the values.
               muImg(itr1, itr2) = aggMu;
               varImg(itr1, itr2) = aggSigma;
          end
     end
end

