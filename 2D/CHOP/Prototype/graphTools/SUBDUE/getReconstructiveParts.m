%> Name: getReconstructiveParts
%>
%> Description: Given a set of parts in bestSubs, this function greedily
%> selects a set of parts that minimize the likelihood of the data. The data
%> is grouped into overlapping receptive fields, and the reduction in the
%> cost is associated with increasing likelihood of the underlying data. Two
%> factors are contributing towards the data likelihood description, namely
%> node label and position prediction. 
%> 
%> @param bestSubs Initial set of substructures.
%> @param realNodeLabels The real labels of underlying data. 
%> @param realEdgeLabels The real labels of the edges that encode spatial
%> distributions in the bottom level. 
%> @param allEdges All edges encoded in the first level, with each cell
%> corresponding to a separate node's edges. 
%> @param allEdgeProbs Probabilities associated with edges.
%> @param numberOfFinalSubs Selection will stop if the number of selected
%> subs exceeds numberOfFinalSubs. 
%> @param stoppingCoverage The minimum coverage that is required to stop
%> selection.
%> @param uniqueChildren The ids of the nodes(data) to be covered.
%>
%> @retval validSubs Ids of final subs.
%> @retval overallCoverage The coverage on the data.
%> @retval dataLikelihood Data likelihood given the selected parts.
%>
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 08.10.2015
function [validSubs] = getReconstructiveParts(bestSubs, numberOfFinalSubs, level1Coords, uniqueChildren, allLeafNodes)

   coverageStoppingVal = 0.97;
   minLimit = 0.01;
   numberOfBestSubs = numel(bestSubs);
   remainingFirstLevelNodes = unique(cat(2, allLeafNodes{uniqueChildren}));
%   numberOfChildren = {bestSubs.edges};
%   numberOfChildren = cellfun(@(x) size(x,1) + 1, numberOfChildren).^(1/3);
%   numberOfChildren = cellfun(@(x) size(x,1) + 1, numberOfChildren);
%   numberOfChildren = cellfun(@(x) size(x,1) + 1, numberOfChildren);
   firstGraphNodeCount = numel(remainingFirstLevelNodes);
   
   % Find number of leaf nodes per image.
   leafNodeImageIds = single(level1Coords(:,1));
   allImageIds = 1:max(leafNodeImageIds);
   imageLeafNodeCounts = hist(leafNodeImageIds(remainingFirstLevelNodes), allImageIds);
   
   % For images which do not have any leaf nodes left, we assign the count
   % 1 to prevent division by zero.
   imageLeafNodeCounts(imageLeafNodeCounts == 0) = 1;
   numberOfRemainingImages = numel(unique(leafNodeImageIds(remainingFirstLevelNodes)));
   
   % Allocate space for log likelihood results.
   subCoveredNodes = cell(numberOfBestSubs,1);
   
   % Go over all possible part-subpart pairs, and calculate probabilities.
   parfor subItr = 1:numberOfBestSubs
       instanceChildren = bestSubs(subItr).instanceChildren;
       
        % Save child probabilities.
        allNodes = fastsortedunique(sort(instanceChildren));
        allNodes = (allNodes(:))';
        coveredNodes = fastsortedunique(sort(cat(2, allLeafNodes{allNodes})));
        subCoveredNodes{subItr} = coveredNodes;
   end
   childrenCounts = cellfun(@(x) numel(x), subCoveredNodes);
   childrenCountLimits = minLimit * childrenCounts;
   
     %% Finally, we implement an algorithm for part selection. 
     % This step is done to perform an initial pass to reduce the number of
     % parameters (subs to be selected) significiantly. Then, we perform
     % another pass using a data likelihood measure. This step implements a
     % coverage-based part selection mechanism.
     subCounter = 0; 
     addedValueArr = [];
     selectedSubIdx = zeros(1,numberOfBestSubs) > 0;
     markedNodes = zeros(size(level1Coords,1),1) > 0;
     valueArr = inf(1,numberOfBestSubs);
     
     % We calculate image-based coverage and then round it.
     validNodeArr = ones(size(level1Coords,1),1) > 0;

     % Mark invalid nodes.
     invalidArr = cellfun(@(x) numel(x), subCoveredNodes);
     valueArr(invalidArr == 0) = 0;
     subLimit = numberOfFinalSubs;
     cumVal = 0;
     displayThr = 0;
     while subCounter < subLimit
       maxLocVal = -inf;
       maxLoc = 0;
       maxSubIdx = [];
       maxLocMarkedNodes = [];

       % Choose eligible subs to choose from, and loop over them.
       eligibleSubs = find(valueArr);
       curSubItr = 1;
       while curSubItr <= numel(eligibleSubs)
           subItr = eligibleSubs(curSubItr);
           tempMarkedNodes = markedNodes;
           tempSubIdx = selectedSubIdx;
           tempSubIdx(subItr) = 1;
           children = subCoveredNodes{subItr};
           
           % Remove already defined children!
           children = children(~markedNodes(children) & validNodeArr(children));
           subCoveredNodes{subItr} = children;

           % Calculate value of a sub.
           tempMarkedNodes(children) = 1;
 %          diffVal = (nnz(tempMarkedNodes) - prevVal) * numberOfChildren(subItr);
 %          diffVal = (nnz(tempMarkedNodes) - prevVal) / childrenCounts(subItr);
 %          diffVal = ((nnz(tempMarkedNodes) - prevVal) * numberOfChildren(subItr))/ childrenCounts(subItr);
           diffVal = numel(children);
 
           % Save diffVal.
           if diffVal < valueArr(subItr)
                if diffVal < childrenCountLimits(subItr)
                     valueArr(subItr) = 0;
      %               valueArr(subItr) = diffVal;
                else
                     valueArr(subItr) = diffVal;
                end
           end

           % Save value if this part has maximum value.
           if diffVal > 0 && diffVal > maxLocVal
                if maxLocVal > -inf && subCounter > 0
                    curSubItr = 1;
                    eligibleSubs = find(valueArr > maxLocVal);
                end
                maxLocVal = diffVal;
                maxLocMarkedNodes = tempMarkedNodes;
                maxLoc = subItr;
                maxSubIdx = tempSubIdx;
           else
                curSubItr = curSubItr + 1;
           end
      end
      if isempty(maxSubIdx)
           break;
      end

      % Save info, and move on to the next iteration.
      addedValueArr = [addedValueArr, maxLocVal]; %#ok<AGROW>
      valueArr(maxSubIdx) = 0;
      markedNodes = maxLocMarkedNodes;
      cumVal = cumVal + maxLocVal;
      selectedSubIdx = maxSubIdx;
      subCounter = subCounter + 1;
      
      % Calculate coverage. 
      imageCoverages = hist(leafNodeImageIds(markedNodes), allImageIds);
      achievedImages = imageCoverages ./ imageLeafNodeCounts >= coverageStoppingVal;
      validNodeArr(ismembc(leafNodeImageIds, single(find(achievedImages)))) = 0;
      
      % Calculate coverage, and check if we've covered enough data.
      % Then, break if necessary.
%       coverage = nnz(markedNodes) / firstGraphNodeCount;
%        if coverage >= coverageStoppingVal 
%            break;
%        end
      if nnz(achievedImages) == numberOfRemainingImages
         break; 
      end

      % Print output.
      if nnz(achievedImages)/numberOfRemainingImages >= displayThr
            display(['Selected  sub # ' num2str(subCounter) ' with id ' ...
               num2str(maxLoc) ', and ' num2str(nnz(achievedImages)) ' out of ' num2str(numberOfRemainingImages) ' images have been covered at least %' num2str(coverageStoppingVal*100) '.']);
         displayThr = displayThr + 0.1;
      end
     end
        
        
    validSubs = find(selectedSubIdx >= 0.5);
    subCoveredNodes = subCoveredNodes(validSubs);
   
   % Calculate statistics.
   allCoveredNodes = cat(2, subCoveredNodes{:});
   allCoveredNodes = unique(allCoveredNodes);
   overallCoverage = numel(allCoveredNodes) / firstGraphNodeCount;
   
   % Printing.
   display(['[SUBDUE] We have selected  ' num2str(numel(validSubs)) ...
        ' out of ' num2str(numel(bestSubs)) ' subs.. Coverage: ' num2str(overallCoverage) '.']);

end
