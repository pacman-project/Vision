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

    numberOfBestSubs = numel(bestSubs);
    allChildren = cell(numberOfBestSubs,1);
    for bestSubItr = 1:numel(bestSubs)
        instanceValidationIdx = bestSubs(bestSubItr).instanceValidationIdx > 0;
        children = bestSubs(bestSubItr).instanceChildren(instanceValidationIdx, :);
        
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
    allChildren = unique(cat(1, allChildren{:}));
    prevGraphNodeCount = numel(allChildren);

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
   
   display(['[SUBDUE] Starting threshold search with ' num2str(midThr) '. We are limited to ' num2str(hiddenNodeCount) ' compositions.']);
   display(['[SUBDUE] Initially, we have ' num2str(numberOfBestSubs) ' subs to consider.']);
   %% Searching for optimal threshold using a binary search mechanism.
   while (currentDepth <= maxDepth)
       %% Optimality logic is implemented here. We're searching for an optimal threshold.
       % We have a flag to indicate if we have close-to-optimal coverage. 
       isSolutionOptimal = false;

       % Calculate which leaf nodes are covered.
       finalSubList = zeros(numberOfBestSubs, 1);
       subNodes = cell(numberOfBestSubs,1);
       subMatchScores = cell(numberOfBestSubs,1);
       for bestSubItr = 1:numberOfBestSubs
           adaptiveThreshold = ((midThr * (size(bestSubs(bestSubItr).edges,1) * 2 + 1)) + singlePrecision);
           validInstanceIdx = bestSubs(bestSubItr).instanceMatchCosts < adaptiveThreshold;
           validInstanceIdx = bestSubs(bestSubItr).instanceValidationIdx > 0 & validInstanceIdx;
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
       validSubIdx = getDisjointSubs(bestSubs, ...
           nodeDistanceMatrix, edgeDistanceMatrix, singlePrecision, midThr);
       finalSubList(validSubIdx) = 1;
       numberOfRemainingBestSubs = nnz(validSubIdx);
       
       % Eliminating invalid subs by setting their contributions to zero.
       % We select all subs, while maximizing our metric.
       remainingNodes = cat(1, subNodes{validSubIdx});
       if ~isempty(remainingNodes)
            remainingNodes = fastsortedunique(sort(remainingNodes));
       else
            remainingNodes = [];
       end
       
       % We calculate our optimality metric. For now, it's unsupervised,
       % and gives equal weight to coverage/mean match scores.
       overallCoverage = numel(remainingNodes) / prevGraphNodeCount;
       overallMatchScore = sum(cellfun(@(x) sum(x), subMatchScores(validSubIdx)));
       overallMatchScoreDenom = sum(cellfun(@(x) numel(x), subMatchScores(validSubIdx)));
       if overallMatchScoreDenom ~= 0
            overallMatchScore = overallMatchScore / overallMatchScoreDenom;
       end
      
       subMeanCoverages = cellfun(@(x) numel(x) / prevGraphNodeCount, subNodes);
       subMeanMatchScores = cellfun(@(x) mean(x), subMatchScores);
       
       % This is the multiplication to find our metric.
       unsupMetric = overallCoverage * overallMatchScore;
       
       % Printing.
       display(['[SUBDUE] We have selected  ' num2str(numberOfRemainingBestSubs) ...
            ' out of ' num2str(numberOfBestSubs) ' subs.. Coverage * average match cost: ' num2str(unsupMetric) '.']);
       
       % If the score is ideal, we mark solution flag as true.
       if unsupMetric > stoppingCoverage
           isSolutionOptimal = true;
       end
       
       % Update final part list.
       finalSubList = finalSubList(finalSubList > 0);
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
            elseif ~isSolutionOptimal && abs(overallCoverage - overallMatchScore) <= singlePrecision
                display(['[SUBDUE] Found an OK similarity threshold (best we can do): ' num2str(midThr) '! Quitting..']);
                break;
            elseif ~isSolutionOptimal && overallCoverage < overallMatchScore
                    minThr = midThr;
                    midThr = (maxThr + minThr) / 2;
                    display(['[SUBDUE] Not optimal coverage. ' ...
                    ' Increasing the threshold to ' num2str(midThr) ' and continuing to search..']);
            elseif ~isSolutionOptimal && overallCoverage > overallMatchScore
                maxThr = midThr;
                midThr = (maxThr + minThr) / 2;
                display(['[SUBDUE] Matching cost too high.' ...
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