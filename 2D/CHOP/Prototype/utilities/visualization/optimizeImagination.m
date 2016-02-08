%> Name: optimizePoE
%>
%> Description: We project high level nodes to an image, calculate
%> gradients for nodes across all layers, and slowly move them around to
%> climb the data likelihood function. 
%>
%> @param nodes High-level nodes in format of (realLabelId, posX, posY,
%levelItr; ...]. 
%> @param vocabulary Part vocabulary.
%> @param imageSize Image size in pixels.
%> @param options Program options.
%> 
%> @retval optimizedImg Optimized image (backprojected).
%>
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 07.02.2016
function [ optimizedImg ] = optimizeImagination( nodes, vocabulary, imageSize, options )
     stepSize = 1;
     filters = options.filters;
     filters = cellfun(@(x) (x - min(min(x))) / (max(max(x)) - min(min(x))), filters, 'UniformOutput', false);
     
     %% First, imagine layer 1 nodes.
     [experts, parseTrees] = projectNode(nodes, vocabulary, 1, 'modal');
     
    % Center the nodes.
    experts = double(experts);
    minX = min(experts(:,2));
    maxX = max(experts(:,2));
    minY = min(experts(:,3));
    maxY = max(experts(:,3));
    midPoint = round([(minX+maxX)/2, (minY+maxY)/2]);
    experts(:,2:3) = round((experts(:,2:3) - repmat(midPoint, size(experts,1), 1))+ repmat(imageSize/2, size(experts,1), 1));
     
     %% Second, we start a gradient-descent procedure to move things around so they match.
     numberOfNodes = max(max(parseTrees));
     parseTreesSize = size(parseTrees);
     pixelPredictions = cell(imageSize);
     modalImg = zeros(imageSize);
     
     gradients = zeros(parseTreesSize(1),1);
     for nodeItr = 1:numberOfNodes
          idx = find(parseTrees == nodeItr);
          [rows, cols] = ind2sub(parseTreesSize, idx);
          
          % If this is the only node, move on.
          if numel(idx) == size(parseTrees,1)
               continue;
          end
          
          % Get previous likelihood value to compare.
          % Obtain fixed experts that will not move.
          fixedIdx = ones(parseTreesSize(1),1) > 0;
          fixedIdx(rows) = 0;
          fixedExperts = experts(fixedIdx,:);
          
          % Find likelihood of the modal reconstruction.
  %        refLikelihoodVal = findLikelihood(experts, [], imageSize, options, filters);
          [curLikelihoodVal, modalImg, curPixelPredictions] = findLikelihood(fixedExperts, pixelPredictions, modalImg, imageSize, options, filters);
          
          % Calculate gradients by placing this node's leaves in four
          % different directions (up, down, right, left) if possible.
          newExperts = experts(rows,:);
          
          for directionItr = 1:4
               switch directionItr
                    case 1
                         offsets = [-stepSize,0];
                    case 2
                         offsets = [stepSize,0];
                    case 3
                         offsets = [0,-stepSize];
                    case 4
                         offsets = [0,stepSize];
               end
               
               % Update expert positions.
               newExperts(:,2:3) = newExperts(:,2:3) + repmat(offsets, numel(rows),1);
               
               % Calculate gradient.
               
          end
     end
end

function [likelihoodVal, pixelPredictions, modalImg] = findLikelihood(experts, pixelPredictions, modalImg, imgSize, options, filters)
     % Program arguments.
     halfFilterSize = round((size(options.filters{1},1) - 1) / 2);
     maxInfluenceRadius = ceil(halfFilterSize * sqrt(2));
     
     % Eliminate unnecessary nodes.
     experts = experts(experts(:,2) > 1 & experts(:,2) < imgSize(1) + 1 & ...
           experts(:,3) > 1 & experts(:,3) < imgSize(2) + 1, :);

     % Create overlap array.
     overlapMatrix = zeros(imgSize) > 0;
     idx = sub2ind(imgSize, experts(:,2), experts(:,3));
     overlapMatrix(idx) = 1;
     maxOverlapMatrix = imdilate(overlapMatrix, strel('disk', maxInfluenceRadius));
     
     % Find sigma values for every distance.
     maxDist = ceil(sqrt(imgSize(1)^2 + imgSize(2)^2)) + 1;
     distVals = 0:1/(maxInfluenceRadius):1;
     gaussVals = normpdf(distVals, 0, 0.4);
     allSigmaVals = 1./gaussVals;
     allSigmaVals = allSigmaVals(1:(maxInfluenceRadius+1))';
     
     % Calculate remaining program variables.
     posArr = double(experts(:,2:3));
     gaborIdArr = double(experts(:,1));
     numberOfGabors = size(experts,1);
     if size(experts,2) > 3
          activationArr = experts(:,4);
     else
          activationArr = ones(numberOfGabors,1);
     end
     
     firstHalf = ceil(options.gaborFilterSize/2) - 1;
     secHalf = options.gaborFilterSize - (firstHalf+1);
     
     % Calculate pixel-level likelihoods.
     for itr1 = 1:imgSize(1)
          for itr2 = 1:imgSize(2)
               if ~maxOverlapMatrix(itr1, itr2)
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
                              predictionArr(expItr,1) = 0;
                         end
                         predictionArr(expItr,2) = sigmaVals(expId);
                    end
               end
               predictionArr = [pixelPredictions{itr1,itr2}; predictionArr]; %#ok<AGROW>
               
               % Calculate predictions using gaussian mixtures.
               probArr = zeros(256,size(predictionArr,1));
               numberOfExperts = size(predictionArr,1);
               weights = (predictionArr(:,2)') .^ -1;
               weights = weights/sum(weights);
               for gaborItr = 1:numberOfExperts
                    likelihoodArr = mvncdf([(valItr/256) - 1/512; (valItr/256)+1/512], muVals(gaborItr), sigmaVals(gaborItr));
                    probArr(valItr+1,gaborItr) = likelihoodArr(2) - likelihoodArr(1);
               end
%                        probSums = sum(probArr,1);
%                        probArr = probArr ./ repmat(probSums, 256, 1);
               probArr = sum(probArr .* repmat(weights, 256, 1),2);

               % Get the most likely explanation, and
               % calculate its log likelihood.
               [probability, aggMu] = max(probArr);
               likelihoodImg(itr1, itr2) = probability;
               meanImg(itr1, itr2) = aggMu/256;
          end
     end

end
