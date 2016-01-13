function [ muImg, varImg ] = obtainPoE( level1Nodes, imgSize, options )
     imgSize = double(imgSize);
     
     % Find sigma values for every distance.
     maxDist = ceil(sqrt(imgSize(1)^2 + imgSize(2)^2));
     distVals = 0:1/(maxDist+1):1;
     gaussVals = normpdf(distVals, 0, 0.4);
     allSigmaVals = 1./gaussVals;
     allSigmaVals = allSigmaVals(1:(maxDist+1))';
     
     % Calculate remaining program variables.
     posArr = double(level1Nodes(:,2:3));
     gaborIdArr = double(level1Nodes(:,1));
     numberOfGabors = size(level1Nodes,1);
     filters = options.filters;
     filters = cellfun(@(x) (x - min(min(x))) / (max(max(x)) - min(min(x))), filters, 'UniformOutput', false);
     firstHalf = ceil(options.gaborFilterSize/2) - 1;
     secHalf = options.gaborFilterSize - (firstHalf+1);
     
     % Calculate pixel-level likelihoods.
     muImg = zeros(imgSize);
     varImg = zeros(imgSize);

     for itr1 = 1:imgSize(1)
          for itr2 = 1:imgSize(2)
               location = [itr1, itr2];
               distances = repmat(location, numberOfGabors, 1) - posArr;
               actualDistances = round(sqrt(sum(distances.^2,2))) + 1;
               sigmaVals = allSigmaVals(actualDistances);
               overlappingIds = distances(:,1) >(-firstHalf-1) & distances(:,1) < (secHalf+1) & distances(:,2) >(-firstHalf-1) & distances(:,2) < (secHalf+1);

               % If there's at least one overlapping filter, get
               % prediction by multiplying distributions.
               aggSigma = sigmaVals(1);
               if overlappingIds(1) == 1
                    filterId = gaborIdArr(1);
                    filterVals = filters{filterId};
                    validPos = distances(1,:) + [firstHalf, secHalf] + 1;
                    aggMu = filterVals(validPos(1), validPos(2));
               else
                    aggMu = 0;
               end

               for prodItr = 2:numberOfGabors
                    % Get the mu of the new distribution.
                    if overlappingIds(prodItr) == 1
                         filterId = gaborIdArr(prodItr);
                         filterVals = filters{filterId};
                         validPos = distances(prodItr,:) + [firstHalf, secHalf] + 1;
                         newMu = filterVals(validPos(1), validPos(2));   
                    else
                         newMu = 0;
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
end

