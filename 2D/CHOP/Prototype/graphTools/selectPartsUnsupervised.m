%> Name: selectPartsUnsupervised
%>
%> Description: Given a set of parts in bestSubs, this function 
%> calculates the optimal set of parts on the validation data, if it
%> exists. If it does not exist, training data is used. Each part is
%> evaluated based on its coverage on the data, and total matching cost of
%> the part's instances. We're trying to maximize the coverage on the
%> training set, while minimizing the matching cost of the instances. Only
%> a limited set of subs is selected. 
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
%> Ver 1.0 on 13.04.2015
function [bestSubs, optimalThreshold] = selectPartsUnsupervised(bestSubs, ...
    nodeDistanceMatrix, edgeDistanceMatrix, ...
    singlePrecision, stoppingCoverage, numberOfFinalSubs, minThreshold, maxThreshold, maxDepth)

    fitValidationData = 0;
    numberOfBestSubs = numel(bestSubs);
    allChildren = cell(numberOfBestSubs,1);
    for bestSubItr = 1:numel(bestSubs)
        if fitValidationData
            instanceValidationIdx = bestSubs(bestSubItr).instanceValidationIdx > 0;
            children = bestSubs(bestSubItr).instanceChildren(instanceValidationIdx, :);
        else
            children = bestSubs(bestSubItr).instanceChildren;
        end
        
        % Get unique children.
        if ~isempty(children)
              children  = children(:);
              children = fastsortedunique(sort(children));
        else 
              children = [];
        end
        
        % Make the sizes uniform, i.e. Nx1.
        if size(children,2) ~= 1
            children = children';
        end
        
        % Assign the children.
        allChildren{bestSubItr} = children;
    end
    uniqueChildren = fastsortedunique(sort(cat(1, allChildren{:})));
    prevGraphNodeCount = numel(uniqueChildren);

    % Learn the maximum children count in instances.
    maxSubSize = 1;
    for bestSubItr = 1:numberOfBestSubs
        maxSubSize = max(maxSubSize, (size(bestSubs(bestSubItr).edges,1)+1));
    end
        
   % Search variables.
   hiddenNodeCount = numberOfFinalSubs; 
   minThr = minThreshold; % Minimum threshold for elastic matching.
   maxThr = maxThreshold; % Max threshold for elastic part matching. 
   midThr = (minThr + maxThr) / 2;
   currentDepth = 1;
   smartSubElimination = 1;
   
   display(['[SUBDUE] Starting threshold search with ' num2str(midThr) '. We are limited to ' num2str(hiddenNodeCount) ' compositions.']);
   display(['[SUBDUE] Initially, we have ' num2str(numberOfBestSubs) ' subs to consider.']);
   %% Searching for optimal threshold using a binary search mechanism.
   while (currentDepth <= maxDepth)
       %% Optimality logic is implemented here. We're searching for an optimal threshold.
       % We have a flag to indicate if we have close-to-optimal coverage. 
       isSolutionOptimal = false;
       isCoverageOptimal = false;
       
       % Calculate which leaf nodes are covered.
       subNodes = cell(numberOfBestSubs,1);
       subMatchScores = cell(numberOfBestSubs,1);
       for bestSubItr = 1:numberOfBestSubs
           adaptiveThreshold = ((midThr * (size(bestSubs(bestSubItr).edges,1) * 2 + 1)) + singlePrecision);
           validInstanceIdx = bestSubs(bestSubItr).instanceMatchCosts < adaptiveThreshold;
           if fitValidationData
               validInstanceIdx = bestSubs(bestSubItr).instanceValidationIdx > 0 & validInstanceIdx;
           end
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
           validSubs = find(validSubIdx);
           display(['[SUBDUE] Considering only disjoint examples, we are down to ' num2str(numel(validSubs)) ' subs to consider.']);
       else
           validSubs = (1:numel(bestSubs))';
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
       numberOfRemainingBestSubs = numel(validSubs);
       while addedSubOffset <= numberOfFinalSubs
           numberOfCoveredNodes = nnz(nodeArr);
           
           % Print current progress.
           if rem(addedSubOffset,10) == 0
                display(['[SUBDUE] Selecting sub ' num2str(addedSubOffset) '/' num2str(numberOfRemainingBestSubs) ...
                    ' out of ' num2str(numberOfBestSubs) '.. Current coverage: ' num2str(numberOfCoveredNodes/prevGraphNodeCount) '.']);
           end
           
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
            if addedSubOffset > numberOfFinalSubs && ...
                    currentDepth == maxDepth
                numberOfFinalSubs = numberOfFinalSubs + 1;
            end
       end
       validSubs = find(selectedSubIdx);
       
       %% Calculate optimality of the solution.
       % Mark remaining instance nodes.
       numberOfRemainingBestSubs = numel(validSubs);
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

       if overallCoverage >= stoppingCoverage
          isCoverageOptimal = true; 
       end
       
       % This is the multiplication to find our metric.
       unsupMetric = overallCoverage * overallMatchScore;
       
       % Printing.
       display(['[SUBDUE] We have selected  ' num2str(numberOfRemainingBestSubs) ...
            ' out of ' num2str(numberOfBestSubs) ' subs.. Coverage: ' num2str(overallCoverage) ', average match score:' num2str(overallMatchScore) '.']);
       
       % If the score is ideal, we mark solution flag as true.
       if unsupMetric > stoppingCoverage
           isSolutionOptimal = true;
       end
       
       % Update final part list.
       finalSubList = validSubs;
       optimalThreshold = midThr;
    
       %% Depending on the feedback, we're continuing or stopping the search.
        % Depending on which side we are on the ideal hidden node count, 
        % we make a binary search on the "good" threshold.
        currentDepth = currentDepth + 1;
        if (currentDepth <= maxDepth)
            if isSolutionOptimal
                % If we've found the perfect number of hidden nodes, exit.
                display(['[SUBDUE] Found a perfect similarity threshold: ' num2str(midThr) '! Quitting..']);
                break;
            elseif numel(finalSubList) < numberOfFinalSubs && ~isCoverageOptimal
                if smartSubElimination
                    smartSubElimination = 0;
                    currentDepth = currentDepth - 1;
                    display('[SUBDUE] Not optimal coverage. Disabling smart sub elimination and trying again.');
                else
                    minThr = midThr;
                    midThr = (maxThr + minThr) / 2;
                    display(['[SUBDUE] Not optimal coverage. ' ...
                    ' Increasing the threshold to ' num2str(midThr) ' and continuing to search..']);
                end
            elseif ~isCoverageOptimal
                    minThr = midThr;
                    midThr = (maxThr + minThr) / 2;
                    display(['[SUBDUE] Not optimal coverage. ' ...
                    ' Increasing the threshold to ' num2str(midThr) ' and continuing to search..']);
            else
                maxThr = midThr;
                midThr = (maxThr + minThr) / 2;
                display(['[SUBDUE] Trying to reduce matching costs.' ...
                    ' Lowering the threshold to ' num2str(midThr) ' and continuing to search..']);
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
        sub.instanceValidationIdx = sub.instanceValidationIdx(validInstances,:);
        bestSubs(bestSubItr) = sub;
   end
end