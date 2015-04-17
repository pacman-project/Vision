function [validSubs, overallCoverage, overallMatchScore] = getReconstructiveParts(bestSubs, ...
    numberOfFinalSubs, moreSubsAllowed, smartSubElimination, midThr, stoppingCoverage, uniqueChildren, nodeDistanceMatrix, ...
    edgeDistanceMatrix, singlePrecision)

   numberOfBestSubs = numel(bestSubs);
   prevGraphNodeCount = numel(uniqueChildren);
   
   % Calculate which leaf nodes are covered.
   subNodes = cell(numberOfBestSubs,1);
   subMatchScores = cell(numberOfBestSubs,1);
   for bestSubItr = 1:numberOfBestSubs
       adaptiveThreshold = ((midThr * (size(bestSubs(bestSubItr).edges,1) * 2 + 1)) + singlePrecision);
       validInstanceIdx = bestSubs(bestSubItr).instanceMatchCosts < adaptiveThreshold;
       nodes = bestSubs(bestSubItr).instanceChildren(validInstanceIdx, :);
       nodes = nodes(nodes > 0);
       nodes = nodes(:);
       if ~isempty(nodes)
           nodes = fastsortedunique(sort(nodes));
       else
           nodes = [];
       end
       subNodes(bestSubItr) = {nodes};

       % We calculate match scores of the instances.
       matchCosts = bestSubs(bestSubItr).instanceMatchCosts(validInstanceIdx, :);
       subMatchScores(bestSubItr) = {(adaptiveThreshold - matchCosts) / adaptiveThreshold};
   end

   % Eliminate subs that match to better subs. We'll have subs that are
   % far away from each other (in terms of pairwise distances), and
   % thus we will have less subs to evaluate.]
   if smartSubElimination
       validSubIdx = getDisjointSubs(bestSubs, ...
           nodeDistanceMatrix, edgeDistanceMatrix, singlePrecision, midThr);
%       display(['[SUBDUE] Considering only disjoint examples, we are down to ' num2str(nnz(validSubIdx)) ' subs to consider.']);
   else
       validSubIdx = ones(numel(bestSubs),1) > 0;
   end
   
   %% We select a number of subs that cover most, if not all of the data.
   % Now, we apply a greedy approach here that selects subs one by one,
   % by their contributions. Ideally, we can replace this step with
   % an optimization problem solver such as a genetic algorithm, or
   % simulated annealing.
   maxContributions = cellfun(@(x) numel(x), subNodes);
   maxContributions(~validSubIdx) = 0;
   subMeanMatchScores = cellfun(@(x) mean(x), subMatchScores);
   subMeanMatchScores(isnan(subMeanMatchScores)) = 0;
   nodeArr = zeros(max(uniqueChildren),1) > 0;
   selectedSubIdx = zeros(numberOfBestSubs,1) > 0;
   addedSubOffset = 1;
%   numberOfRemainingBestSubs = numel(validSubs);
   while addedSubOffset <= numberOfFinalSubs
       numberOfCoveredNodes = nnz(nodeArr);

%        % Print current progress.
%        if rem(addedSubOffset,10) == 0
%             display(['[SUBDUE] Selecting sub ' num2str(addedSubOffset) '/' num2str(numberOfRemainingBestSubs) ...
%                 ' out of ' num2str(numberOfBestSubs) '.. Current coverage: ' num2str(numberOfCoveredNodes/prevGraphNodeCount) '.']);
%        end

       % If we've covered enough nodes, stop.
       if numberOfCoveredNodes/prevGraphNodeCount >= stoppingCoverage
           break;
       end

       % Select next best sub.
       valueArr = inf(numberOfBestSubs,1);
       curValue = 0;
       for bestSubItr = 1:numberOfBestSubs
            % If we've already marked this sub or it is not promising, go on.
            if maxContributions(bestSubItr) <= curValue
                continue; 
            end

            % Calculate value of this node.
            tempFlagArr = nodeArr;
            tempFlagArr(subNodes{bestSubItr}) = 1;
            tempValue = nnz(tempFlagArr) - numberOfCoveredNodes;
            tempValue = tempValue * subMeanMatchScores(bestSubItr);

            % Record the value for the end of iteration.
            valueArr(bestSubItr) = tempValue;
            if tempValue > curValue
                curValue = tempValue;
            end
        end
        maxContributions = min(valueArr, maxContributions);
        valueArr(isinf(valueArr)) = 0;

        % If there is a new part that introduces novelty, add the best one
        % to the list of selected subs, and move on.
        [value, maxLoc] = max(valueArr);
        if value > 0
            selectedSubIdx(maxLoc) = 1;
            nodeArr(subNodes{maxLoc}) = 1;
            maxContributions(maxLoc) = 0;
        else
            % No new info can be introduced by any subs, just stop.
            break;
        end
        addedSubOffset = addedSubOffset + 1;

        % If we're at the end of the search, and optimal coverage has
        % not been met, we increase final number of subs.
        if addedSubOffset > numberOfFinalSubs && moreSubsAllowed
            numberOfFinalSubs = numberOfFinalSubs + 1;
        end
   end
   % Record preserved subs.
   validSubs = find(selectedSubIdx);
   
   % Mark remaining instance nodes.
   remainingNodes = cat(1, subNodes{validSubs});
   if ~isempty(remainingNodes)
        remainingNodes = fastsortedunique(sort(remainingNodes));
   else
        remainingNodes = [];
   end
   
   % We calculate our optimality metric. For now, it's unsupervised,
   % and gives equal weight to coverage/mean match scores.
   overallCoverage = numel(remainingNodes) / prevGraphNodeCount;
   overallMatchScore = sum(cellfun(@(x) sum(x), subMatchScores(validSubs)));
   overallMatchScoreDenom = sum(cellfun(@(x) numel(x), subMatchScores(validSubs)));
   if overallMatchScoreDenom ~= 0
        overallMatchScore = overallMatchScore / overallMatchScoreDenom;
   end
   
   % Printing.
%   display(['[SUBDUE] We have selected  ' num2str(numel(validSubs)) ...
%        ' out of ' num2str(numberOfBestSubs) ' subs.. Coverage: ' num2str(overallCoverage) ', average match score:' num2str(overallMatchScore) '.']);

end