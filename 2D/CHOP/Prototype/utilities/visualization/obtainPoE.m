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
function [modalImg, varMat, likelihoodMat, likelihoodVal] = obtainPoE(experts, modalImg, varMat, likelihoodMat, imgSize, filters, likelihoodLookupTable, expertTypes)
     % Program arguments.
     % If added expert flags haven't been given, we consider all experts as new.
     if nargin == 7
          expertTypes = ones(size(experts,1),1);
          prevExpertsGiven = false;
     else
          if nnz(expertTypes == 0) > 0
               prevExpertsGiven = true;
          else
               prevExpertsGiven = false;
          end
     end
     
     % Save which experts are to be processed.
     if prevExpertsGiven
          newExperts = experts(expertTypes == 1,:);
          removedExperts = experts(expertTypes == -1,:);
     else
          newExperts = experts;
          removedExperts = [];
     end
     
     % Calculate program arguments and allocate space for images.
     minPixelValue = 1/255;
     if isempty(modalImg)
          modalImg =  ones(imgSize) * minPixelValue;
          varMat = zeros(imgSize);
          likelihoodMat = zeros(imgSize);
     end
     filterSize = size(filters,1);
     filterMatrixSize = size(filters);
     minSigma = 0.15;
     
     % If fast flag is active, we don't take non-overlapping filters into
     % account.
     maxInfluenceRadius = ceil(filterSize/2);
     maxDist = maxInfluenceRadius^2;
     
     %% Now, we obtain product of non-influencing experts from a range of 1 to max.
     % Eliminate unnecessary nodes.
     numberOfNewExperts = nnz(expertTypes == 1);
     numberOfOldExperts = nnz(expertTypes == 0);
     numberOfRemovedExperts = nnz(expertTypes == -1);
     numberOfFinalExperts = numberOfOldExperts + numberOfNewExperts - numberOfRemovedExperts;
     finalExpertIdx = expertTypes > -1;
     aggSigma = minSigma / sqrt(numberOfFinalExperts);
     
     % Find modified experts to focus on them. 
     changedExperts = [newExperts; removedExperts];
     changedExperts = changedExperts(changedExperts(:,2) > 1 & changedExperts(:,2) < imgSize(1) + 1 & ...
           changedExperts(:,3) > 1 & changedExperts(:,3) < imgSize(2) + 1, :);
     
     % Create overlap array.
     maxOverlapMatrix = zeros(imgSize) > 0;
     idx = sub2ind(imgSize, round(changedExperts(:,2)), round(changedExperts(:,3)));
     maxOverlapMatrix(idx) = 1;
     cx=maxInfluenceRadius+1;cy=cx;ix=2*maxInfluenceRadius+1;iy=ix;r=maxInfluenceRadius;
     [x,y]=meshgrid(-(cx-1):(ix-cx),-(cy-1):(iy-cy));
     c_mask=((x.^2+y.^2)<r^2);
     maxOverlapMatrix = imdilate(maxOverlapMatrix, c_mask);
     
     % We create a max overlap matrix based on the types of the filters.
     % Calculate remaining program variables.
     posArr = double(experts(finalExpertIdx, 2:3));
     gaborIdArr = double(experts(finalExpertIdx, 1));
     
     % Create program variables, and find changed pixels.
     firstHalf = ceil(filterSize/2) - 1;
     secHalf = filterSize - (firstHalf+1);
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
     if numberOfNewExperts ~= numberOfRemovedExperts
          likelihoodMat = (likelihoodMat * numberOfOldExperts + ...
               dummyLog * (numberOfNewExperts - numberOfRemovedExperts)) / numberOfFinalExperts;
     else
          likelihoodMat(:) = dummyLog;
     end
     
     % Calculate pixel-level likelihoods.
     for pixelItr = 1:numberOfPixels
          % Obtain x and y coordinates.
          itr1 = posX(pixelItr);
          itr2 = posY(pixelItr);
          
          % Initialize structures.
          logProb = 0;
          location = [itr1, itr2];
          
          % Finally, we get the filters that overlap with this point.
          relativeLocations = -posArr;
          relativeLocations(:,1) = relativeLocations(:,1) + location(1);
          relativeLocations(:,2) = relativeLocations(:,2) + location(2);
          overlappingIdx = sum(relativeLocations.^2,2) < maxDist;
          overlappingLocations = relativeLocations(overlappingIdx,:);
          overlappingLocations(:,1) = firstHalf + 1 + overlappingLocations(:,1); 
          overlappingLocations(:,2) = secHalf + 1 + overlappingLocations(:,2); 
          
          % Keep track of overlapping gabors filters, and use them in
          % calculations.
          overlappingIdx = find(overlappingIdx);
          overlappingIdx = overlappingIdx(expertTypes(overlappingIdx) > -1);
          overlappingIdxCount = numel(overlappingIdx);

          % No experts? Move on. We're confident.
           predictionIdx = overlappingLocations(:,1) + (overlappingLocations(:,2)-1)*filterMatrixSize(1) + (gaborIdArr(overlappingIdx)-1)*filterMatrixSize(1)*filterMatrixSize(2);
           predictionArr = filters(predictionIdx);
          
          %% Calculate product of experts! We're using equal sigmas, which greatly reduces computation.
          aggMu = round(sum(predictionArr)/overlappingIdxCount);
          
          % Update Mu and Sigma of the final predictions.
          modalImg(itr1, itr2) = aggMu;
          varMat(itr1, itr2) = aggSigma;
                    
          %% If there are values in the prediction array, we move forward.
          % The product calculation gives us normalized values. Now,
          % let's switch back to unnormalized values.
          logProb = logProb + sum(likelihoodLookupTable(predictionArr(:,1) + 1, (aggMu+1)));
          
          % Finally, add non-overlapping predictions.
          noncontributingExpertCount = numberOfFinalExperts - overlappingIdxCount;
          if noncontributingExpertCount > 0
               logProb = logProb + dummyLikelihoodVals(noncontributingExpertCount);
          end

          %% Save the information.
          likelihoodMat(itr1, itr2) = logProb / numberOfFinalExperts;
     end
     
     % Return the likelihood value.
     likelihoodVal = sum(sum(likelihoodMat));
end