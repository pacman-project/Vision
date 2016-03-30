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
function [modalImg, likelihoodMat, likelihoodVal] = obtainPoE(experts, modalImg, likelihoodMat, imgSize, filters, likelihoodLookupTable, expertTypes)
     % Before we do anything, let's find the experts on the edge of the
     % image and eliminate them for fast processing.
     filterSize = size(filters,1);
     maxInfluenceRadius = ceil(filterSize/2);
     validExpertIdx = experts(:,2) >= maxInfluenceRadius & experts(:,3) >= maxInfluenceRadius & ...
         experts(:,2) <= (imgSize(1)-(maxInfluenceRadius-1)) & experts(:,3) < (imgSize(2)-(maxInfluenceRadius-1));
     experts = experts(validExpertIdx,:);
     IA = (1:size(experts,1))';
%     [experts, IA, ~] = unique(experts, 'rows', 'stable');

     % Program arguments.
     % If added expert flags haven't been given, we consider all experts as new.
     if nargin == 6
          expertTypes = ones(size(experts,1),1);
          prevExpertsGiven = false;
     else
          expertTypes = expertTypes(validExpertIdx);
          expertTypes = expertTypes(IA);
          if nnz(expertTypes == 0) > 0
               prevExpertsGiven = true;
          else
               prevExpertsGiven = false;
          end
     end
     likelihoodFlag = ~isempty(likelihoodLookupTable);
     
     % Save which experts are to be processed.
     if prevExpertsGiven
          newExperts = experts(expertTypes == 1,:);
          removedExperts = experts(expertTypes == -1,:);
     else
          newExperts = experts;
          removedExperts = [];
     end
     
     % Calculate program arguments and allocate space for images.
     if isempty(modalImg)
          modalImg =  ones(imgSize, 'uint8');
          if likelihoodFlag
               likelihoodMat = zeros(imgSize);
          end
     end
     filterMatrixSize = size(filters);
     minSigma = 0.15;
     
     %% Now, we obtain product of non-influencing experts from a range of 1 to max.
     % Eliminate unnecessary nodes.
     numberOfNewExperts = nnz(expertTypes == 1);
     numberOfOldExperts = nnz(expertTypes < 1);
     numberOfRemovedExperts = nnz(expertTypes == -1);
     numberOfFinalExperts = numberOfOldExperts - numberOfRemovedExperts + numberOfNewExperts;
     finalExpertIdx = expertTypes > -1;
     
     % Find modified experts to focus on them. 
     changedExperts = [newExperts; removedExperts];
     
     % Create overlap array.
     emptyImageMatrix = zeros(imgSize);
     maxOverlapMatrix = emptyImageMatrix > 0;
     idx = sub2ind(imgSize, round(changedExperts(:,2)), round(changedExperts(:,3)));
     maxOverlapMatrix(idx) = 1;
     cx=maxInfluenceRadius;cy=cx;ix=2*maxInfluenceRadius-1;iy=ix;r=maxInfluenceRadius;
     [x,y]=meshgrid(-(cx-1):(ix-cx),-(cy-1):(iy-cy));
     c_mask=((x.^2+y.^2)<r^2);
     maxOverlapMatrix = imdilate(maxOverlapMatrix, c_mask);
     
     % We create a max overlap matrix based on the types of the filters.
     % Calculate remaining program variables.
     posArr = single(experts(finalExpertIdx, 2:3));
     gaborIdArr = single(experts(finalExpertIdx, 1));
     
     % Create program variables, and find changed pixels.
     changedIdx = find(maxOverlapMatrix);
     [posX, posY] = ind2sub(imgSize, changedIdx);
     numberOfPixels = numel(changedIdx);
     
     %% Create dummy likelihood values for non-contributing experts.
     dummyProb = normcdf([-1/510, 1/510], 0, minSigma);
     dummyProb = dummyProb(2) - dummyProb(1);
     dummyLog = log(dummyProb);
     dummyLikelihoodVals = (1:numberOfFinalExperts) * dummyLog;
     
     % Create dummy likelihood and variance values and add them to the
     % computation.
     if likelihoodFlag && numberOfNewExperts ~= numberOfRemovedExperts
          likelihoodMat = (likelihoodMat * numberOfOldExperts + ...
               dummyLog * (numberOfNewExperts - numberOfRemovedExperts)) / numberOfFinalExperts;
     elseif likelihoodFlag
          likelihoodMat(:) = dummyLog;
     end
     
     distances = pdist2(posArr,single([posX, posY])) < maxInfluenceRadius;
     allExpertIds = 1:numberOfFinalExperts;
     negPosArr = -posArr;
     
     % Calculate pixel-level likelihoods.
     for pixelItr = 1:numberOfPixels
          % Obtain x and y coordinates.
          itr1 = posX(pixelItr);
          itr2 = posY(pixelItr);
          
          % Initialize structures.
          logProb = 0;
          location = [itr1, itr2] + maxInfluenceRadius;
          logIdx = distances(:, pixelItr);
          overlappingIdx = allExpertIds(logIdx);
           
          % Keep track of overlapping gabors filters, and use them in
          % calculations.
          overlappingIdxCount = numel(overlappingIdx);
          
          if overlappingIdxCount > 0
               % Obtain relative locations.
               overlappingLocations = negPosArr(logIdx,:);
               overlappingLocations(:,1) = overlappingLocations(:,1) + location(1);
               overlappingLocations(:,2) = overlappingLocations(:,2) + location(2);

               % No experts? Move on. We're confident.
               predictionIdx = overlappingLocations(:,1) + (overlappingLocations(:,2)-1)*filterMatrixSize(1) + (gaborIdArr(overlappingIdx)-1)*filterMatrixSize(1)*filterMatrixSize(2);
               predictionArr = filters(predictionIdx);
               
               %% Calculate product of experts! We're using equal sigmas, which greatly reduces computation.
               aggMu = round(sum(predictionArr)/overlappingIdxCount);

               % Calculate data likelihood.
               if likelihoodFlag
                    logProb = logProb + sum(likelihoodLookupTable(predictionArr(:,1) + 1, (aggMu+1)));
               end
          else
               aggMu = 1;
          end
          
          % Update Mu and Sigma of the final predictions.
          modalImg(itr1, itr2) = uint8(aggMu);
         
          % Calculate likelihoods if needed.
          if likelihoodFlag
               % Finally, add non-overlapping predictions.
               noncontributingExpertCount = numberOfFinalExperts - overlappingIdxCount;
               if noncontributingExpertCount > 0
                    logProb = logProb + dummyLikelihoodVals(noncontributingExpertCount);
               end
               
               %% Save the information.
               likelihoodMat(itr1, itr2) = logProb / numberOfFinalExperts;
          end
     end
     
     % Return the likelihood value.
     if likelihoodFlag
          likelihoodVal = sum(sum(likelihoodMat));
     else
          likelihoodVal = 0;
     end
end