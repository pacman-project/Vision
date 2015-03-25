%> Name: getReconstructiveParts
%>
%> Description: Given a set of parts in bestSubs, it finds a minimal
%> subset that can represent the training data with least redundancy.
%> Normally, this should be considered as an optimization problem operating
%> in very large space (bestSubs can have tens of thousands of parts).
%> However, since we do not have a smart algorithm that can guide us in this
%> huge search space, we opt for a "sort-of" depth-first methodology. It looks 
%> like this:
%> 
%> Given a set of potential parts in bestSubs = {P_1, ...,  P_N}
%> 1) Select P_i which covers most non-covered area on training set
%> 2) Mark the areas exhibited by P_i as 'covered',
%> 3) Exclude P_i from bestSubs
%> 4) Go back to step 1 until %stoppingCoverage of the data is covered, or no new 
%> areas (leaf nodes) can be covered. 
%>
%> @param bestSubs A set of substructures evaluated on the training set.
%> @param instanceChildrenDescriptors N_instance x maxNumberOfChildren
%> array that has an ordered list of children from the previous level for
%> each instance.
%> @param remainingInstanceLabels N_instance x 1 array that marks the part
%> label for each instance.
%> @param allLeafNodes N x 1 cell array that keeps what leaf nodes are covered 
%> by each part's instances.  
%> prevGraphNodeCount Number of nodes in the previous level's object
%> graphs.
%> @param stoppingCoverage Selection will stop if stoppingCoverage percent
%> of previous levels' leaf nodes are already covered.
%> @param numberOfFinalSubs Selection will stop if the number of selected
%> subs exceeds numberOfFinalSubs. 
%>
%> @retval bestSubs Final set of best substructures.
%> @retval instanceChildrenDescriptors N_instance x maxNumberOfChildren
%> array that has an ordered list of children from the previous level for
%> each instance.
%> @retval remainingInstanceLabels N_instance x 1 array that marks the part
%> label for each instance.
%>
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 23.02.2015
function [bestSubs, optimalThreshold] = getReconstructiveParts(bestSubs, ...
    allLeafNodes, prevGraphNodeCount, maxPrevGraphNodeId, nodeDistanceMatrix, edgeDistanceMatrix, ...
    singlePrecision, stoppingCoverage, numberOfFinalSubs, minThreshold, maxThreshold, maxDepth)

    numberOfBestSubs = numel(bestSubs);
    % Learn the maximum children count in instances.
    maxSubSize = 1;
    for bestSubItr = 1:numberOfBestSubs
        maxSubSize = max(maxSubSize, (size(bestSubs(bestSubItr).edges,1)+1));
    end
        
   % Search variables.
   hiddenNodeCount = numberOfFinalSubs; 
   minThr = minThreshold; % Minimum threshold for elastic matching.
   maxThr = maxThreshold; % Max threshold for elastic part matching. 
   midThr = (minThr + maxThr) / 2 + singlePrecision;
   currentDepth = 1;
   
   % Put leaf nodes in a different form.
   maxLeafNodeCount = max(cellfun(@(x) numel(x), allLeafNodes));
   leafNodeArr = zeros(numel(allLeafNodes), maxLeafNodeCount, 'int32');
   for nodeItr = 1:numel(allLeafNodes)
       leafNodeRow = allLeafNodes{nodeItr};
       leafNodeArr(nodeItr, 1:numel(leafNodeRow)) = leafNodeRow;
   end
   
   display(['[SUBDUE] Starting threshold search with ' num2str(midThr) '. We are limited to ' num2str(hiddenNodeCount) ' compositions.']);
   display(['[SUBDUE] Initially, we have ' num2str(numberOfBestSubs) ' subs to consider.']);
   %% Searching for optimal threshold using a binary search mechanism.
   while (currentDepth <= maxDepth)
       %% Optimality logic is implemented here. We're searching for an optimal threshold.
       % We have a flag to indicate if we have close-to-optimal coverage. 
       isCoverageOptimal = false;

        % We have an array to mark whether each leaf node has been detected or not.
       nodeFlagArr = zeros(maxPrevGraphNodeId,1) > 0;
       numberOfBestSubs = numel(bestSubs);

       % Calculate which leaf nodes are covered.
       finalSubList = zeros(numberOfFinalSubs, 1);
       finalSubListNovelty = zeros(numberOfFinalSubs, 1);
       addedSubs = 1;
       subLeafNodes = cell(numberOfBestSubs,1);
       for bestSubItr = 1:numberOfBestSubs
           validInstanceIdx = bestSubs(bestSubItr).instanceMatchCosts < ((midThr * (size(bestSubs(bestSubItr).edges,1) * 2 + 1)) + singlePrecision);
           allLeafNodes = leafNodeArr(bestSubs(bestSubItr).instanceChildren(validInstanceIdx, :), :);
           allLeafNodes = allLeafNodes(:);
           subLeafNodes(bestSubItr) = {allLeafNodes(allLeafNodes>0)};
       end
       
       % partLeafCounts holds the latest info for all parts. It tracks the maximum 
       % number of novel leaf nodes can introduce to the system by adding it, 
       % but the actual number of leaf nodes that are new is likely to be lower.
       partLeafCounts = cellfun(@(x) numel(x), subLeafNodes); 
       [partLeafCounts, partLeafCountOrder] = sort(partLeafCounts, 'descend');
       curSubLeafNodes = subLeafNodes(partLeafCountOrder);

       % Eliminate subs that match to better subs. We'll have subs that are
       % far away from each other (in terms of pairwise distances), and
       % thus we will have less subs to evaluate.]
       validSubIdx = getDisjointSubs(bestSubs(partLeafCountOrder), ...
           nodeDistanceMatrix, edgeDistanceMatrix, singlePrecision, midThr);
       numberOfRemainingBestSubs = nnz(validSubIdx);
%       numberOfRemainingBestSubs = numberOfBestSubs;
       display(['[SUBDUE] Eliminating overlapping subs. We are down to ' num2str(numberOfRemainingBestSubs) ' subs.']);
       % Eliminating invalid subs by setting their contributions to zero.
       partLeafCounts(~validSubIdx) = 0;
       
       % Go ahead and select parts.
       for bestSubItr = 1:numberOfFinalSubs
            % Get the contribution of this sub in terms of number of leaf
            % nodes.
            tempFlagCount = nnz(nodeFlagArr);
            if rem(bestSubItr,10) == 0
                display(['[SUBDUE] Selecting sub ' num2str(bestSubItr) '/' num2str(numberOfRemainingBestSubs) ...
                    ' out of ' num2str(numberOfBestSubs) '.. Current coverage: ' num2str(tempFlagCount/prevGraphNodeCount) '.']);
            end

            % The stopping criterion is set to covering stoppingCoverage percent of all available
            % leaf nodes.
            if tempFlagCount/prevGraphNodeCount >= stoppingCoverage            
               isCoverageOptimal = true;
               break; 
            end

            % Compare the contribution of this node to the rest, and pick
            % the next best one.
            valueArr = inf(numberOfBestSubs,1);
            novelLeafNodes = zeros(numberOfBestSubs,1);
            curValue = 0;
            for bestSubItr2 = 1:numberOfBestSubs
                % If we've already marked this sub or it is not promising, go on.
                if partLeafCounts(bestSubItr2) <= curValue
                    continue; 
                end

                % Calculate value of this node.
                tempFlagArr = nodeFlagArr;
                tempFlagArr(curSubLeafNodes{bestSubItr2}) = 1;
                tempValue = nnz(tempFlagArr) - tempFlagCount;

                % Record the value for the end of iteration.
                novelLeafNodes(bestSubItr2) = tempValue;
                valueArr(bestSubItr2) = tempValue;
                curValue = tempValue;
            end
            partLeafCounts = min(valueArr, partLeafCounts);
            valueArr(isinf(valueArr)) = 0;

            % If there is a new part that introduces novelty, add the best one
            % to the list of selected subs, and move on.
            [value, maxLoc] = max(valueArr);
            if value > 0
                finalSubList(addedSubs) = maxLoc;
                numberOfNovelLeafNodes = novelLeafNodes(maxLoc);
                finalSubListNovelty(bestSubItr) = numberOfNovelLeafNodes / numel(curSubLeafNodes{maxLoc});
                nodeFlagArr(curSubLeafNodes{maxLoc}) = 1;
                partLeafCounts(maxLoc) = 0;
                addedSubs = addedSubs+1;
            else
                % No new info can be introduced by any subs, just stop.
                isCoverageOptimal = true;
                break;
            end
       end
       
       
       display(['[SUBDUE] We have selected  ' num2str(bestSubItr) ...
            ' out of ' num2str(numberOfBestSubs) ' subs.. Coverage: ' num2str(nnz(nodeFlagArr)/prevGraphNodeCount) '.']);
       
       % If everything is covered, we mark coverage flag as true.
       if nnz(nodeFlagArr) == prevGraphNodeCount
           isCoverageOptimal = true;
       end
       
       % Update final part list.
       finalSubList = finalSubList(finalSubList > 0);
       finalSubList = sort(finalSubList);
       finalSubList = partLeafCountOrder(finalSubList);
       optimalThreshold = midThr;
    
       %% Depending on the feedback, we're continuing or stopping the search.
        % Depending on which side we are on the ideal hidden node count, 
        % we make a binary search on the "good" threshold.
        currentDepth = currentDepth + 1;
        if (currentDepth <= maxDepth)
            if numel(finalSubList) == hiddenNodeCount && isCoverageOptimal
                % If we've found the perfect number of hidden nodes, exit.
                display(['[SUBDUE] Found perfect similarity threshold: ' num2str(midThr) '! Quitting..']);
                break;
            elseif numel(finalSubList) < hiddenNodeCount
                maxThr = midThr;
                midThr = (maxThr + minThr) / 2;
                display(['[SUBDUE] Too few generated compositions (' num2str(numel(finalSubList)) ').' ...
                    ' Lowering the threshold to ' num2str(midThr) ' and continuing to search..']);
            else
                % We have few hidden nodes. Make threshold
                % larger.                              
                minThr = midThr;
                midThr = (maxThr + minThr) / 2;
                display(['[SUBDUE] Too many generated compositions (' num2str(numel(finalSubList)) ').' ...
                    ' Increasing the threshold to ' num2str(midThr) ' and continuing to search..']);
            end
        end
   end

   % Update instance information.
   bestSubs = bestSubs(finalSubList);
   % Update bestSubs instances by taking the new threshold into account.
   for bestSubItr = 1:numel(bestSubs)
        sub = bestSubs(bestSubItr);
        validInstances = sub.instanceMatchCosts < ((size(sub.edges,1)*2+1) * optimalThreshold + singlePrecision);
        sub.instanceCenterIdx = sub.instanceCenterIdx(validInstances,:);
        sub.instanceChildren = sub.instanceChildren(validInstances,:);
        if ~isempty(sub.edges)
            sub.instanceEdges = sub.instanceEdges(validInstances,:);
        end
        sub.instanceSigns = sub.instanceSigns(validInstances,:);
        sub.instanceCategories = sub.instanceCategories(validInstances,:);
        sub.instanceMatchCosts = sub.instanceMatchCosts(validInstances,:);
        bestSubs(bestSubItr) = sub;
   end
end