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
    singlePrecision, stoppingCoverage, numberOfFinalSubs, minThreshold, ...
    maxThreshold, maxDepth, validationFolds)
    
    % Keep best subs.
    orgBestSubs = bestSubs;
    
    % Start processing with the remaining subs.
    numberOfOrgBestSubs = numel(bestSubs);
    allChildren = cell(numberOfOrgBestSubs,1);
    for bestSubItr = 1:numel(bestSubs)
        children = bestSubs(bestSubItr).instanceChildren;

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
    
    %% We're performing cross-validation. 
    % Each time, we're excluding subs that have been enumerated from a
    % validation subset, and using the rest to cover the data. After we
    % find the thresholds for eachs set, we simply take their average to
    % find an optimal threshold.
    crossValThresholds = zeros(validationFolds,1, 'single');
    subEliminationFlags = zeros(validationFolds,1);
    parfor valItr = 1:validationFolds
        % We exclude the subs which have zero-cost matchs on this subset,
        % but not on other subsets.
        validSubIdx = ones(numel(orgBestSubs),1) > 0;
        
        % If there's only one fold, or we do not wish to do
        % cross-validation, we fit the training data.
        if validationFolds > 1
            for subItr = 1:numel(orgBestSubs)
                instanceMatchCosts = orgBestSubs(subItr).instanceMatchCosts;
                instanceValidationIdx = orgBestSubs(subItr).instanceValidationIdx;

                if nnz(instanceMatchCosts == 0 & instanceValidationIdx ~= valItr) == 0
                   validSubIdx(subItr) = 0; 
                end
            end
        end
        bestSubs = orgBestSubs(validSubIdx);
        numberOfBestSubs = numel(bestSubs);
        uniqueChildren = fastsortedunique(sort(cat(1, allChildren{validSubIdx}))); %#ok<PFBNS>

        % Learn the maximum children count in instances.
        maxSubSize = 1;
        for bestSubItr = 1:numberOfBestSubs
            maxSubSize = max(maxSubSize, (size(bestSubs(bestSubItr).edges,1)+1));
        end

       % Search variables.
       minThr = minThreshold; % Minimum threshold for elastic matching.
       maxThr = maxThreshold; % Max threshold for elastic part matching. 
       midThr = (minThr + maxThr) / 2;
       optimalThreshold = midThr;
       currentDepth = 1;
       moreSubsAllowed = 0;
       smartSubElimination = 1;

       display(['[SUBDUE] Starting threshold search with ' num2str(midThr) '. We are limited to ' num2str(numberOfFinalSubs) ' compositions.']);
       display(['[SUBDUE] Initially, we have ' num2str(numberOfBestSubs) ' subs to consider.']);
       %% Searching for optimal threshold using a binary search mechanism.
       while (currentDepth <= maxDepth)
           % If this is the final leg of search, we allow part selection to
           % have more subs than the max number.
           if currentDepth == maxDepth
              moreSubsAllowed = 1;
           end
           
           %% Optimality logic is implemented here. We're searching for an optimal threshold.
           % First, we obtain non-overlapping, minimal set of subs.
           [validSubs, isSolutionOptimal, isCoverageOptimal, overallCoverage, overallMatchScore] = getReconstructiveParts(bestSubs, ...
                numberOfFinalSubs, moreSubsAllowed, smartSubElimination, midThr, ...
                stoppingCoverage, uniqueChildren, nodeDistanceMatrix, ...
                edgeDistanceMatrix, singlePrecision);
           
           % If the coverage is good, mark the corresponding flag.
           if overallCoverage >= stoppingCoverage
              isCoverageOptimal = true; 
           end

           % This is the multiplication to find our metric.
           unsupMetric = overallCoverage * overallMatchScore;

           % If the score is ideal, we mark solution flag as true.
           if unsupMetric > stoppingCoverage
               isSolutionOptimal = true;
           end

           % Update final part list.
           optimalThreshold = midThr;

           %% Depending on the feedback, we're continuing or stopping the search.
            % Depending on which side we are on the ideal hidden node count, 
            % we make a binary search on the "good" threshold.
            currentDepth = currentDepth + 1;
            if (currentDepth <= maxDepth)
                if isSolutionOptimal
                    % If we've found the perfect number of hidden nodes, exit.
                    display(['[SUBDUE] Found a perfect similarity threshold: ' num2str(midThr) '! Quitting..']);
                elseif numel(validSubs) < numberOfFinalSubs && ~isCoverageOptimal
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
       crossValThresholds(valItr) = optimalThreshold;
       subEliminationFlags(valItr) = smartSubElimination;
    end
    
    %% We have found the optimal threshold. 
    % Now, we obtain the subs one final time using the new threshold, and exit.
    bestSubs = orgBestSubs;
    if validationFolds > 1
        display(['[SUBDUE] Thresholds learned from cross-validation folds: ' ...
            mat2str(crossValThresholds) ', with mean: ' num2str(mean(crossValThresholds)) '.']);
    end
    optimalThreshold = mean(crossValThresholds);
    smartSubElimination = round(mean(subEliminationFlags));
    uniqueChildren = fastsortedunique(sort(cat(1, allChildren{:})));
    
    % Get a new set of subs.
   [finalSubList, ~, ~, ~, ~] = getReconstructiveParts(bestSubs, ...
    numberOfFinalSubs, 1, smartSubElimination, optimalThreshold, ...
    stoppingCoverage, uniqueChildren, nodeDistanceMatrix, ...
    edgeDistanceMatrix, singlePrecision);

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