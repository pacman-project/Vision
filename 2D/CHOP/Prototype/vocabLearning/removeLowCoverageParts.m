function [allCliques, allInstancePositions, newMetaData, newChildrenArr, newPartDistributions] = removeLowCoverageParts(allCliques, allInstancePositions, newMetaData, newChildrenArr, newPartDistributions, prevLevel, leafNodes, level1Coords, nodeCoverageThr, halfRFSize)
     subsetCount = 100;
    
     % Eliminate parts which have low coverage. 
     [~, firstInstances, ~] = unique(newMetaData(:,1), 'R2012a');
     numberOfParts = numel(firstInstances);
     firstInstances = cat(1, firstInstances, size(newMetaData,1)+1);
     
     % Decide on image ids for every clique.
     imageIds = prevLevel(allCliques(:,end) , 4);
     
      % We wish to learn different feature maps here.
      % Start by calculating instance centers.
      allLeafCovers = cell(numberOfParts,1);
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
            tmpCount = min(subsetCount, ceil(size(relevantCliques, 1)/2));
            selectedCliques = datasample(1:size(relevantCliques, 1), tmpCount);
            
            % Calculate covered nodes.
            for instanceItr = selectedCliques
                 tmpArr = cat(2, leafNodes{relevantCliques(instanceItr,:)});
                 tmpNodes = fastsortedunique(sort(tmpArr));
                 leafNodeCounts(instanceItr) = numel(tmpNodes);
 %                leafCovers{instanceItr} = tmpNodes;
            end
  %          allLeafCovers{partItr} = leafCovers;
            allLeafCounts{partItr} = leafNodeCounts;
      end
%     allLeafCovers = cat(1, allLeafCovers{:}); %#ok<NASGU>
     realNodeCovers = cat(1, allLeafCounts{:});
       
     % Calculate max coverage factor for all instances here.
     comparedIdx = realNodeCovers > 0;
     discretePositions = round(allInstancePositions(comparedIdx, :));
     imageIds = imageIds(comparedIdx);
     realNodeCovers = realNodeCovers(comparedIdx, :);
     processedImages = fastsortedunique(sort(imageIds));
     [~, imageIdx, ~] = unique(level1Coords(:,1), 'R2012a');
     imageIdx = cat(1, imageIdx, size(level1Coords,1)+1);
     maxNodeCovers = zeros(size(discretePositions,1), 1, 'single');
     
     for imageItr = 1:numel(processedImages)
          tmpIdx = imageIds == processedImages(imageItr);
          relevantPositions = level1Coords(imageIdx(processedImages(imageItr)):(imageIdx(processedImages(imageItr)+1)-1), 2:3);
          comparedPositions = discretePositions(tmpIdx, :);
          
          % Calculate # of layer 1 nodes in this receptive field. 
          distances = pdist2(relevantPositions, comparedPositions, 'chebychev') < halfRFSize;
          nodeCovers = sum(distances, 1);
          maxNodeCovers(tmpIdx) = nodeCovers;
     end
     
     % Calculate coverage percentages.
     nodeCoverPct = realNodeCovers./maxNodeCovers;
     [~, firstInstances, ~] = unique(newMetaData(comparedIdx,1), 'R2012a');
     firstInstances = cat(1, firstInstances, nnz(comparedIdx)+1);
     
     % Decide which parts to eliminate here.
     partCoverPct = zeros(numberOfParts,1, 'single');
     for partItr = 1:numberOfParts
          partCoverPct(partItr) = mean(nodeCoverPct(firstInstances(partItr):(firstInstances(partItr+1)-1)));
     end
     validParts = find(partCoverPct >= nodeCoverageThr);
     validInstances = ismembc(newMetaData(:,1), int32(validParts));
     
     % Finally, eliminate parts which have low coverage on the overall. 
     allCliques = allCliques(validInstances, :);
     newMetaData = newMetaData(validInstances, :);
     allInstancePositions = allInstancePositions(validInstances, :);
     newChildrenArr = newChildrenArr(validParts, :);
     newPartDistributions = newPartDistributions(validParts, :);
     [~, ~, newMetaData(:,1)] = unique(newMetaData(:,1));
end

