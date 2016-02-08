function [ muImg, varImg ] = calculateGradients( level1Nodes, imgSize, options, fastFlag, muImgInitial, varImgInitial )
     imgSize = double(imgSize);
     filterSize = size(options.filters{1},1);
     halfFilterSize = round((size(options.filters{1},1) - 1) / 2);
     maxInfluenceRadius = halfFilterSize + 1;
     
     % If initial mu/var images are provided, we're only supposed to update
     % what we have.
     if nargin == 4
          muImg = zeros(imgSize);
          varImg = zeros(imgSize);
     else
          muImg = muImgInitial;
          varImg = varImgInitial;
     end
     
     if fastFlag
          % Eliminate unnecessary nodes.
          level1Nodes = level1Nodes(level1Nodes(:,2) > 1 & level1Nodes(:,2) < imgSize(1) + 1 & ...
                level1Nodes(:,3) > 1 & level1Nodes(:,3) < imgSize(2) + 1, :);

          % Create overlap array.
          overlapMatrix = zeros(imgSize) > 0;
          idx = sub2ind(imgSize, level1Nodes(:,2), level1Nodes(:,3));
          overlapMatrix(idx) = 1;
          maxOverlapMatrix = imdilate(overlapMatrix, strel('disk', maxInfluenceRadius));
          overlapMatrix = imdilate(overlapMatrix, strel('square', filterSize));
     else
          overlapMatrix = ones(imgSize) > 0;
          maxOverlapMatrix = ones(imgSize) > 0;
     end
     
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
               actualDistances = round(sqrt(sum(distances.^2,2))) + 1;
               actualDistances(actualDistances > maxDist) = maxDist;
               invalidIdx = actualDistances > maxInfluenceRadius;
               sigmaVals = zeros(numel(actualDistances),1);
               sigmaVals(~invalidIdx) = allSigmaVals(actualDistances(~invalidIdx));
               overlappingIds = distances(:,1) >(-firstHalf-1) & distances(:,1) < (secHalf+1) & distances(:,2) >(-firstHalf-1) & distances(:,2) < (secHalf+1);

               aggMu = muImg(itr1, itr2);
               aggSigma = max(allSigmaVals);
               if varImg(itr1, itr2) > 0
                    aggSigma = varImg(itr1, itr2);
               end
               
               for prodItr = 1:numberOfGabors
                    if (fastFlag && ~overlappingIds(prodItr)) || invalidIdx(prodItr)
                         continue;
                    end
                    
                    % Get the mu of the new distribution.
                    if overlappingIds(prodItr) == 1
                         filterId = gaborIdArr(prodItr);
                         filterVals = filters{filterId};
                         validPos = distances(prodItr,:) + [firstHalf, secHalf] + 1;
                         newMu = filterVals(validPos(1), validPos(2)) * activationArr(prodItr);   
                    else
                         newMu = muImg(itr1, itr2);
                    end
                    % Take product of gaussians.
                    aggMu = (aggMu * (sigmaVals(prodItr)^2) + newMu * (aggSigma^2)) / (aggSigma^2 + sigmaVals(prodItr)^2);
                    aggSigma = sqrt(((aggSigma^2)*(sigmaVals(prodItr)^2))/(aggSigma^2+sigmaVals(prodItr)^2));
               end

               % Assign the values.
               muImg(itr1, itr2) = aggMu;
               varImg(itr1, itr2) = aggSigma;
          end
     end
%     muImg = muImg/max(max(muImg));
end

