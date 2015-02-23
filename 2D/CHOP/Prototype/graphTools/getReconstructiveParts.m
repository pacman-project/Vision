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
%> @param bestSubs Final set of best substructures.
%> @param instanceChildrenDescriptors N_instance x maxNumberOfChildren
%> array that has an ordered list of children from the previous level for
%> each instance.
%> @param remainingInstanceLabels N_instance x 1 array that marks the part
%> label for each instance.
%>
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 23.02.2015
function [bestSubs, instanceChildrenDescriptors, remainingInstanceLabels] = getReconstructiveParts(bestSubs, ...
    instanceChildrenDescriptors, remainingInstanceLabels, allLeafNodes, prevGraphNodeCount, stoppingCoverage, numberOfFinalSubs)

    % We have an array to mark whether each leaf node has been detected or not.
   nodeFlagArr = zeros(prevGraphNodeCount,1) > 0;
   numberOfBestSubs = numel(bestSubs);

   % Calculate which leaf nodes are covered.
   finalSubList = zeros(numberOfFinalSubs, 1);
   addedSubs = 1;
   subLeafNodes = cell(numberOfBestSubs,1);
   for bestSubItr = 1:numberOfBestSubs
       subChildren = setdiff(unique(instanceChildrenDescriptors(remainingInstanceLabels == bestSubItr, :)),0);
       subLeafNodes(bestSubItr) = {unique([allLeafNodes{subChildren}])};
   end
   
   % partLeafCounts holds the latest info for all parts. It tracks the maximum 
   % number of novel leaf nodes can introduce to the system by adding it, 
   % but the actual number of leaf nodes that are new is likely to be lower.
   partLeafCounts = cellfun(@(x) numel(x), subLeafNodes); 

   % Go ahead and select parts.
   for bestSubItr = 1:numberOfFinalSubs
        % Get the contribution of this sub in terms of number of leaf
        % nodes.
        tempFlagCount = nnz(nodeFlagArr);
        if rem(bestSubItr,10) == 0
            display(['[SUBDUE] Selecting sub ' num2str(bestSubItr) '/' num2str(numberOfFinalSubs) ...
                ' out of ' num2str(numberOfBestSubs) '.. Current coverage: ' num2str(tempFlagCount/numel(nodeFlagArr)) '.']);
        end

        % The stopping criterion is set to covering stoppingCoverage percent of all available
        % leaf nodes.
        if tempFlagCount/numel(nodeFlagArr) >= stoppingCoverage
           break; 
        end
        
        % Compare the contribution of this node to the rest, and pick
        % the next best one.
        valueArr = inf(numberOfBestSubs,1);
        curValue = 0;
        for bestSubItr2 = 1:numberOfBestSubs
            % If we've already marked this sub or it is not promising, go on.
            if partLeafCounts(bestSubItr2) <= curValue
                continue; 
            end
           
            % Calculate value of this node.
            tempFlagArr = nodeFlagArr;
            tempFlagArr(subLeafNodes{bestSubItr2}) = 1;
            tempValue = nnz(tempFlagArr) - tempFlagCount;

            % Record the value for the end of iteration.
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
            nodeFlagArr(subLeafNodes{maxLoc}) = 1;
            partLeafCounts(maxLoc) = 0;
            addedSubs = addedSubs+1;
        else
            % No new info can be introduced by any subs, just stop.
            break;
        end
   end
   % Update final part list.
   finalSubList = finalSubList(finalSubList > 0);
   finalSubList = sort(finalSubList);

   % Update instance information.
   remainingIdx = ismember(remainingInstanceLabels, finalSubList);
   instanceChildrenDescriptors = instanceChildrenDescriptors(remainingIdx,:);
   remainingInstanceLabels = remainingInstanceLabels(remainingIdx);
   [~, ~, remainingInstanceLabels] = unique(remainingInstanceLabels, 'stable');
   bestSubs = bestSubs(finalSubList);
end