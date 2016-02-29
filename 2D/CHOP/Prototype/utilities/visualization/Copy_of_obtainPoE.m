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
function [modalImg, varMat, likelihoodMat, likelihoodVal] = obtainPoE(experts, modalImg, varMat, likelihoodMat, imgSize, filters, addedExperts)
     % Program arguments.
     % If added expert flags haven't been given, we consider all experts as
     % new.
     if nargin == 6
          addedExperts = ones(size(experts,1),1) > 0;
          prevExpertsGiven = false;
     else
          if nnz(~addedExperts) > 0
               prevExpertsGiven = true;
          else
               prevExpertsGiven = false;
          end
     end
     
     % Save which experts are to be processed.
     if prevExpertsGiven
          newExperts = experts(addedExperts,:);
     else
          newExperts = experts;
     end
     
     % Calculate program arguments and allocate space for images.
     minPixelValue = 1/255;
     if isempty(modalImg)
          modalImg =  ones(imgSize) * minPixelValue;
          varMat = zeros(imgSize);
          likelihoodMat = zeros(imgSize);
     end
     filterSize = size(filters{1},1);
%     maxSigma = 5;
     minSigma = 0.15;
%     meanSigma = 0.2;
     
     % If fast flag is active, we don't take non-overlapping filters into
     % account.
     maxInfluenceRadius = floor(filterSize/2);
     maxDist = maxInfluenceRadius^2;
     
     % Find sigma values for every distance.
     vals = sqrt(0:maxDist);
     vals = vals / floor(filterSize/2);
 %    allSigmaVals = sigmf(vals,[40 0.9]) * maxSigma;
 %    allSigmaVals(allSigmaVals < minSigma) = minSigma;
     allSigmaVals = repmat(minSigma, size(vals,1), size(vals,2));
     
     % Calculate discrete pixel values.
     pixelCenterValues = (0:255)/255;

     %% Now, we obtain product of non-influencing experts from a range of 1 to max.
     % Eliminate unnecessary nodes.
     newExperts = newExperts(newExperts(:,2) > 1 & newExperts(:,2) < imgSize(1) + 1 & ...
           newExperts(:,3) > 1 & newExperts(:,3) < imgSize(2) + 1, :);

     % Create overlap array.
     maxOverlapMatrix = zeros(imgSize) > 0;
     idx = sub2ind(imgSize, round(newExperts(:,2)), round(newExperts(:,3)));
     maxOverlapMatrix(idx) = 1;
     cx=maxInfluenceRadius+1;cy=cx;ix=2*maxInfluenceRadius+1;iy=ix;r=maxInfluenceRadius;
     [x,y]=meshgrid(-(cx-1):(ix-cx),-(cy-1):(iy-cy));
     c_mask=((x.^2+y.^2)<r^2);
     maxOverlapMatrix = imdilate(maxOverlapMatrix, c_mask);
     
     % We create a max overlap matrix based on the types of the filters.
     % Calculate remaining program variables.
     posArr = double(experts(:,2:3));
     gaborIdArr = double(experts(:,1));
     numberOfExperts = size(experts,1);
     if size(experts,2) > 3
          activationArr = experts(:,4);
     else
          activationArr = ones(numberOfExperts,1);
     end
     
     firstHalf = ceil(filterSize/2) - 1;
     secHalf = filterSize - (firstHalf+1);
     changedIdx = find(maxOverlapMatrix);
     [posX, posY] = ind2sub(imgSize, changedIdx);
     numberOfPixels = numel(changedIdx);
     assignedSigmaVals = repmat(allSigmaVals(end), numberOfExperts, 1);
     
     %% Create dummy likelihood values for non-contributing experts.
%     dummyProb = 1/256;
     dummyProb = normcdf([-1/510, 1/510], 0, minSigma);
     dummyProb = dummyProb(2) - dummyProb(1);
     dummyLog = log(dummyProb);
     dummyLikelihoodVals = (1:numberOfExperts) * dummyLog;
     
     % Create dummy likelihood and variance values and add them to the
     % computation.
     likelihoodMat(:) = likelihoodMat(:) + dummyLog * nnz(addedExperts);
     
     % Calculate pixel-level likelihoods.
     for pixelItr = 1:numberOfPixels
          % Obtain x and y coordinates.
          itr1 = posX(pixelItr);
          itr2 = posY(pixelItr);
          
          % Initialize structures.
          logProb = 0;
          location = [itr1, itr2];
          
          % Get relative positions of experts and calculate distances to each expert.
          distances = repmat(location, numberOfExperts, 1) - posArr;
          actualDistances = sum(distances.^2,2);
          actualDistances(actualDistances>=maxDist) = maxDist;
          distances = round(distances);
          
          % Finally, we get the filters that overlap with this point.
          overlappingIds = distances(:,1) >(-firstHalf-1) & distances(:,1) < (secHalf+1) & distances(:,2) >(-firstHalf-1) & distances(:,2) < (secHalf+1);
          
          % Remove far-away experts.
          validIdx = actualDistances < maxDist;
          
          % Obtain variances.
          idx = round(actualDistances) + 1;
          sigmaVals = assignedSigmaVals;
          sigmaVals(validIdx) = allSigmaVals(idx(validIdx));
          
          % Keep track of overlapping gabors filters, and use them in
          % calculations.
          overlappingIdx = find(validIdx);

          % No experts? Move on. We're confident.
          if isempty(overlappingIdx)
               predictionArr = [];
          else 
               predictionArr = zeros(numel(overlappingIdx), 2);

               %% Calculate product of experts.
               for expItr = 1:numel(overlappingIdx)
                    expId = overlappingIdx(expItr);
                    % Get the mu of the new distribution.
                    if overlappingIds(expId) == 1
                         filterId = gaborIdArr(expId);
                         filterVals = filters{filterId};
                         validPos = distances(expId,:) + [firstHalf, secHalf] + 1;
                         predictionArr(expItr,1) = filterVals(validPos(1), validPos(2)) * activationArr(expId);
                    else
                         predictionArr(expItr,1) = minPixelValue;
                    end
                    
                    % Calculate sigma values based on the 2D
                    % gaussians.
                    predictionArr(expItr,2) = sigmaVals(expId);
               end
          end
          
          %% Calculate product of experts.
          % If all experts have the same mu, we can speed up
          % calculation by only calculating sigmas of the multiplied
          % distributions (normalized pdfs).
          if ~isempty(predictionArr)
               uniqueMu = fastsortedunique(sort(predictionArr(:,1)));
               if numel(uniqueMu) == 1
                   aggMu = uniqueMu;                    
                   aggSigma = predictionArr(1,2);
                    for expItr = 2:size(predictionArr,1)
                         newSigma = predictionArr(expItr,2);
                         % Only calculate sigma values.
                         aggSigma = sqrt(((aggSigma^2)*(newSigma^2))/(aggSigma^2+newSigma^2));
                    end
              else
                    aggMu = predictionArr(1,1);
                    aggSigma = predictionArr(1,2);
                    for expItr = 2:size(predictionArr,1)
                         newMu = predictionArr(expItr,1);
                         newSigma = predictionArr(expItr,2);
                         % Take product of gaussians.
                         aggMu = (aggMu * (newSigma^2) + newMu * (aggSigma^2)) / (aggSigma^2 + newSigma^2);
                         aggSigma = sqrt(((aggSigma^2)*(newSigma^2))/(aggSigma^2+newSigma^2));
                    end
               end
               
               % Now, let's convert aggMu to pixel values.
               diffVals = abs(pixelCenterValues - aggMu);
               [~, minIdx] = min(diffVals);
               aggMu = pixelCenterValues(minIdx);
               
               % Update Mu and Sigma of the final predictions.
               modalImg(itr1, itr2) = aggMu;
               varMat(itr1, itr2) = aggSigma;
          end
          
          %% If there are values in the prediction array, we move forward.
          if ~isempty(predictionArr)
               % The product calculation gives us normalized values. Now,
               % let's switch back to unnormalized values.
               vals = repmat([aggMu - 1/510, aggMu + 1/510], size(predictionArr,1), 1);

               %% Calculate log probability of individual distributions at modal point.
               tempVals = normcdf(vals, repmat(predictionArr(:,1),1,2), repmat(predictionArr(:,2),1,2));
               normVals = normcdf(repmat([-1/510, 1+1/510], size(vals,1), 1), repmat(predictionArr(:,1),1,2), repmat(predictionArr(:,2),1,2));
               normVals = normVals(:,2) - normVals(:,1);
               probs = tempVals(:,2) - tempVals(:,1);
               probs = probs ./ normVals;
               logProb = logProb + sum(log(probs));
          end
          
          % Finally, add non-overlapping predictions.
          noncontributingExpertCount = numberOfExperts - numel(overlappingIdx);
          if noncontributingExpertCount > 0
               logProb = logProb + dummyLikelihoodVals(noncontributingExpertCount);
          end

          %% Save the information.
          likelihoodMat(itr1, itr2) = logProb;
     end
     
     likelihoodVal = sum(sum(likelihoodMat));
end


