function [allCliques, allInstancePositions, newMetaData, newChildrenArr, newPartDistributions] = removeLowCoverageParts(allCliques, allInstancePositions, pooledInstancePositions, newMetaData, newChildrenArr, newPartDistributions, prevLevel, leafNodes, level1Coords, nodeCoverageThr, rfSize)
     % Program variables.
     rfCenter = floor(rfSize / 2) + 1;
     lowHalf = 1 - rfCenter;
     highHalf = rfSize - rfCenter;
     maxNodeCount = nnz(level1Coords(:,1) == 1);
     dummyArr = zeros(maxNodeCount,1) > 0;
    
     % Eliminate parts which have low coverage. 
     [~, firstInstances, ~] = unique(newMetaData(:,1), 'R2012a');
     numberOfParts = numel(firstInstances);
     firstInstances = cat(1, firstInstances, size(newMetaData,1)+1);
     
     % Decide on image ids for every clique.
     imageIds = prevLevel(allCliques(:,end) , 4);
     [~, imageOffsets, ~] = unique(level1Coords(:,1), 'R2012a');
     imageOffsets = imageOffsets - 1;
     
     for itr = 1:numel(leafNodes);
          tmp = leafNodes{itr};
          tmp = tmp - imageOffsets(level1Coords(tmp(1)));
          leafNodes{itr} = tmp;
     end
     
      % We wish to learn different feature maps here.
      % Start by calculating instance centers.
      allLeafCounts = cell(numberOfParts,1);
      parfor partItr = 1:numberOfParts
            % Get relevant instances.
            validSubPartCount = nnz(newChildrenArr(partItr,:));
            relevantCliques = allCliques(firstInstances(partItr):(firstInstances(partItr+1)-1),:);
            relevantCliques = relevantCliques(:, (end + (1 - validSubPartCount)):end);
            
            % Calculate leaf covers.
%            leafCovers = cell(size(relevantCliques, 1), 1);
            leafNodeCounts = zeros(size(relevantCliques, 1), 1);
            
            % Randomly select a subset to represent the part.
%            tmpCount = min(subsetCount, ceil(size(relevantCliques, 1)/2));
 %           selectedCliques = datasample(1:size(relevantCliques, 1), tmpCount, 'Replace', false);
            selectedCliques = 1:size(relevantCliques,1);
            % Calculate covered nodes.
            for instanceItr = selectedCliques
                 t1 = dummyArr;
                 t1([leafNodes{relevantCliques(instanceItr,:)}]) = 1;
                 leafNodeCounts(instanceItr) = nnz(t1);
 %                leafCovers{instanceItr} = tmpNodes;
            end
  %          allLeafCovers{partItr} = leafCovers;
            allLeafCounts{partItr} = leafNodeCounts;
      end
%     allLeafCovers = cat(1, allLeafCovers{:}); %#ok<NASGU>
     realNodeCovers = cat(1, allLeafCounts{:});
     
     % Calculate max coverage factor for all instances here.
     processedImages = fastsortedunique(sort(imageIds));
     [~, imageIdx, ~] = unique(level1Coords(:,1), 'R2012a');
     imageIdx = cat(1, imageIdx, size(level1Coords,1)+1);
     maxNodeCovers = cell(numel(processedImages), 1);
     
     % Sort image ids and positions for faster access and parallelization.
     [sortedImageIds, sortedIdx] = sort(imageIds);
     sortedPooledPositions = pooledInstancePositions(sortedIdx,:);
     [~, firstInstances, ~] = unique(sortedImageIds, 'R2012a');
     firstInstances = cat(1, firstInstances, size(sortedImageIds,1)+1);
     
     parfor imageItr = 1:numel(processedImages)
          tmpIdx = firstInstances(imageItr):(firstInstances(imageItr+1)-1);
          relevantPositions = level1Coords(imageIdx(processedImages(imageItr)):(imageIdx(processedImages(imageItr)+1)-1), 4:5);
          comparedPositions = sortedPooledPositions(tmpIdx, :);
          % Check what we have in the receptive field.
          distances = zeros(size(relevantPositions,1), size(comparedPositions,1)) > 0;
          for pointItr = 1:size(comparedPositions,1)
               comparedPoint = comparedPositions(pointItr, :);
               tempPosX = relevantPositions(:, 1) - comparedPoint(1);
               tempPosY = relevantPositions(:, 2) - comparedPoint(2);
               rfCheck = tempPosX >= lowHalf & tempPosX <= highHalf & tempPosY >= lowHalf & tempPosY <= highHalf;
               distances(:, pointItr) = rfCheck;
          end
          
          nodeCovers = sum(distances, 1)';
          maxNodeCovers{imageItr} = nodeCovers;
     end
     maxNodeCovers = cat(1, maxNodeCovers{:});
     maxNodeCovers(sortedIdx) = maxNodeCovers;
     
     % Calculate coverage percentages.
     nodeCoverPct = realNodeCovers./maxNodeCovers;
%     [~, firstInstances, ~] = unique(newMetaData(comparedIdx,1), 'R2012a');
%     firstInstances = cat(1, firstInstances, nnz(comparedIdx)+1);
     
%      % Decide which parts to eliminate here.
%      partCoverPct = zeros(numberOfParts,1, 'single');
%      for partItr = 1:numberOfParts
%           partCoverPct(partItr) = mean(nodeCoverPct(firstInstances(partItr):(firstInstances(partItr+1)-1)));
%      end
%     validParts = find(partCoverPct >= nodeCoverageThr);
%     validInstances = ismembc(newMetaData(:,1), int32(validParts));
     validInstances = nodeCoverPct >= nodeCoverageThr;
     validParts = unique(newMetaData(validInstances,1));
     
     % Finally, eliminate parts which have low coverage on the overall. 
     allCliques = allCliques(validInstances, :);
     newMetaData = newMetaData(validInstances, :);
     allInstancePositions = allInstancePositions(validInstances, :);
     newChildrenArr = newChildrenArr(validParts, :);
     newPartDistributions = newPartDistributions(validParts, :);
     [~, ~, newMetaData(:,1)] = unique(newMetaData(:,1));
end


