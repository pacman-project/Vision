function [validSubs, overallCoverage, overallMatchCost] = getReconstructiveParts(bestSubs, allEdges, allEdgeProbs, ...
    numberOfFinalSubs, valItr, moreSubsAllowed, smartSubElimination, midThr, stoppingCoverage, uniqueChildren, nodeDistanceMatrix, ...
    edgeDistanceMatrix, singlePrecision)

   minNodeProbability = 0.00001;
   numberOfBestSubs = numel(bestSubs);
   prevGraphNodeCount = numel(uniqueChildren);
   maxChildId = max(uniqueChildren);
   prevGraphNodeLogProbs = zeros(1, maxChildId, 'single');
   minLogProb = single(log(minNodeProbability));
   validNodes = zeros(maxChildId, 1) > 0;
   allEdgesArr = {allEdges.adjInfo};
   nonzeroIdx = cellfun(@(x) ~isempty(x), allEdgesArr);
   allEdgesPeripheralNodes = cell(size(allEdgesArr));
   allEdgesPeripheralNodes(nonzeroIdx) = cellfun(@(x) x(:,2), allEdgesArr(nonzeroIdx), 'UniformOutput', false);
   
   %% To start with, we calculate the probability contributions of each sub to the data.
   subLogProbs = cell(numberOfBestSubs, 1);
   subCoveredNodes = cell(numberOfBestSubs,1);
   for subItr = 1:numberOfBestSubs
        % We obtain a partitioning of the data, by getting non-overlapping
        % instances.
       instanceCenterIdx = bestSubs(subItr).instanceCenterIdx;
       [~, IA, ~] = unique(instanceCenterIdx, 'stable');
       numberOfInstances = size(instanceCenterIdx,1);
       validInstances = zeros(numberOfInstances,1) > 0;
       validInstances(IA) = 1;
%       matchedNodes = zeros(maxChildId,1)>0;
%        for instItr = 1:numberOfInstances
%            children = allEdgesPeripheralNodes{instanceCenterIdx(instItr)};
%            if nnz(matchedNodes(children)) == 0
%                matchedNodes(children) = 1;
%                validInstances(instItr) = 1;
%            end
%        end
        
        logProbs = prevGraphNodeLogProbs;
        centerNodes = bestSubs(subItr).instanceCenterIdx(validInstances);
        
        % We consider the center nodes as being reconstructed with great
        % precision.
        coveredNodes = unique(centerNodes);
        logProbs(coveredNodes) = -minLogProb;
        
        if ~isempty(bestSubs(subItr).edges)
             instanceEdges = bestSubs(subItr).instanceEdges(validInstances, :);
             instanceEdges = mat2cell(instanceEdges, ones(size(instanceEdges,1),1), size(instanceEdges,2));
             % Obtain peripheral nodes, their probabilities, and unexplained
             % nodes for this sub. The unexplained nodes are the ones in the RF
             % of this sub's instances, which are not part of the instance
             % definitions.
             peripheralNodes = cellfun(@(x, y) x(y,2), allEdgesArr(centerNodes)', instanceEdges, 'UniformOutput', false);
             peripheralNodes = cat(1, peripheralNodes{:});
             allPeripheralNodes = cat(1, allEdgesArr{centerNodes});
             allPeripheralNodes = allPeripheralNodes(:,2);
             unexplainedNodes = setdiff(allPeripheralNodes, peripheralNodes);
             peripheralNodeProbs = cellfun(@(x, y) x(y), allEdgeProbs(centerNodes)', instanceEdges, 'UniformOutput', false);
             peripheralNodeProbs = cat(1, peripheralNodeProbs{:});

             % Now, we obtain the max probabilities for each node. 
             [vals, idx] = sort(peripheralNodeProbs, 'descend');
             [peripheralNodes, IA, ~] = unique(peripheralNodes(idx), 'stable');
             peripheralNodeLogProbs = -minLogProb + log(vals(IA))';

             % Update log probs.
             logProbs(peripheralNodes) = max(peripheralNodeLogProbs, logProbs(peripheralNodes));
             coveredNodes = unique([coveredNodes; peripheralNodes; unexplainedNodes]);
        end
        logProbs = logProbs(coveredNodes);
        
        % Save the info.
        subCoveredNodes(subItr) = {coveredNodes'};
        subLogProbs(subItr) = {logProbs};
   end
   
   %% Start the main loop. We select subs one by one.
   % This is a greedy selection algorithm. Each new sub's evaluation
   % depends on the previous set, i.e. it shows how much novelty it can bring to the table.
   coveredNodes = zeros(1, maxChildId) > 0;
   selectedSubIdx = zeros(numberOfBestSubs,1) > 0;
   selectedSubs = [];
   curLogProbs = prevGraphNodeLogProbs;
   addedSubOffset = 1;
   reconstructionArr = zeros(1, maxChildId) > 0;
   while addedSubOffset <= numberOfFinalSubs
       % Select next best sub.
       valueArr = -inf(numberOfBestSubs,1);
       curValue = -inf;
       tempCoveredNodes = coveredNodes;
       tempLogProbs = curLogProbs;
       for bestSubItr = 1:numberOfBestSubs
            % If this sub has already been selected, we move on.
            if selectedSubIdx(bestSubItr)
               continue;
            end
            
            % Measure the value of this sub.
            % Get newly covered nodes.
            newCoveredNodes = subCoveredNodes{bestSubItr};
            newCoveredNodesArr = zeros(size(coveredNodes))>0;
            newCoveredNodesArr(newCoveredNodes) = 1;
            newCoveredNodes = newCoveredNodesArr;
            
            % Exclude already covered nodes from new covered nodes.
            validNodeArr = ~reconstructionArr(newCoveredNodes);
            newCoveredNodes = newCoveredNodes & ~reconstructionArr;
            
            % Assign (log) probabilities for these nodes.
            newLogProbs = subLogProbs{bestSubItr};
            newLogProbs = newLogProbs(validNodeArr);
            combinedCoveredNodes = newCoveredNodes | coveredNodes;
            combinedLogProbs = curLogProbs;
            combinedLogProbs(newCoveredNodes) = max(newLogProbs, curLogProbs(newCoveredNodes));
            
            % Calculate value of the sub. We calculate the contribution of
            % this sub to the overall data likelihood, and normalize the
            % value with the square root of the number of covered nodes.
            % This helps single out subs that have very few instances,
            % while providing good coverage on that area.
            contribution = sum(combinedLogProbs(combinedCoveredNodes)) - sum(curLogProbs(combinedCoveredNodes));
%            tempValue = ceil(contribution / sqrt(nnz(newCoveredNodes)));
           tempValue = ceil(contribution);
            
            % Record the value for the end of iteration.
            valueArr(bestSubItr) = tempValue;
            if tempValue > curValue
                maxLoc = bestSubItr;
                tempCoveredNodes = combinedCoveredNodes;
                tempLogProbs =combinedLogProbs;
                curValue = tempValue;
            end
       end
       
       % Mark covered nodes for stoppage criterion.
       prevCoverageCount = nnz(reconstructionArr);
       reconstructionArr(tempLogProbs>0) = 1;
%       reconstructionArr(bestSubs(maxLoc).instanceChildren) = 1;
       
       
       % Update covered nodes and log probs.
       coveredNodes = tempCoveredNodes;
       curLogProbs = tempLogProbs;
       
        % If there is a new part that introduces novelty, add the best one
        % to the list of selected subs, and move on.
        selectedSubIdx(maxLoc) = 1;
        selectedSubs = [selectedSubs, maxLoc];
        addedSubOffset = addedSubOffset + 1;
        
        if rem(addedSubOffset-1,10) == 0 && addedSubOffset > 1
             display(['Selected part ' num2str(maxLoc) ' as the ' num2str(addedSubOffset-1) ...
                  'th part. Its contribution to the overall log likelihood was ' ...
                  num2str(valueArr(maxLoc)) '. It introduced ' num2str(nnz(reconstructionArr) - prevCoverageCount) ...
                  ' new nodes. We are at ' num2str(nnz(reconstructionArr)) ' nodes (out of ' num2str(prevGraphNodeCount) ' possible).']);
        end
        
        % Stopping criterion.
 %        if nnz(curLogProbs > minLogProb) >= stoppingCoverage * prevGraphNodeCount && bestContribution < abs(minLogProb)
         if nnz(reconstructionArr) >= stoppingCoverage * prevGraphNodeCount
              break;
         end
        
        % If we're at the end of the search, and optimal coverage has
        % not been met, we increase final number of subs.
%         if addedSubOffset > numberOfFinalSubs && moreSubsAllowed
%             numberOfFinalSubs = numberOfFinalSubs + 1;
%         end
   end
   overallCoverage = nnz(reconstructionArr) / prevGraphNodeCount;
   overallMatchCost = abs(mean(curLogProbs)/ minLogProb);
   
   % Record preserved subs.
   validSubs = find(selectedSubIdx);
   
   % Printing.
   display(['[SUBDUE] We have selected  ' num2str(numel(validSubs)) ...
        ' out of ' num2str(numberOfBestSubs) ' subs.. Coverage: ' num2str(overallCoverage) ', average normalized match cost:' num2str(overallMatchCost) '.']);

end