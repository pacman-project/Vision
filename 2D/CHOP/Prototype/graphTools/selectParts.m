%> Name: selectParts
%>
%> Description: Given a set of parts in bestSubs, this function 
%> calculates the optimal set of parts on the validation data, if it
%> exists. If it does not exist, training data is used. There are two 
%> evaluation mechanisms for parts: First, we evaluate them using a
%> reconstruction-based criterion. However, if supervisionFlag is true, or
%> our part selection has reduced the categorization performance on the
%> dataset, we switch to discrimination-based part selection.
%> Reconstruction-based selection: Each part is
%> evaluated based on its coverage on the data, and total matching cost of
%> the part's instances. We're trying to maximize the coverage on the
%> training set, while minimizing the matching cost of the instances. Only
%> a limited set of subs is selected. 
%> Discrimination-based selection: We categorize the images, and calculate 
%> the overall accuracy. If the overall accuracy has dropped since the last
%> level, we switch to discrimination, and set supervisionFlag as true.
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
function [bestSubs, optimalThreshold, optimalAccuracy] = selectParts(bestSubs, ...
    nodeDistanceMatrix, edgeDistanceMatrix, ...
    singlePrecision, stoppingCoverage, numberOfFinalSubs, minThreshold, ...
    maxThreshold, maxDepth, validationFolds, validationIdx, categoryArrIdx, imageIdx, supervisionFlag)
    
    % Warn the user. This may take a while.
    display('[SUBDUE] Running unsupervised part selection. This may take a while..');

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
    crossValAccuracy = zeros(validationFolds,1);
    subEliminationFlags = zeros(validationFolds,1);
    if validationFolds > 1
        display(['[SUBDUE/Parallel] Determining the optimal matching threshold with ' num2str(validationFolds) '-fold cross-validation.']); 
    else
        display('[SUBDUE] Finding an optimal matching threshold based on the training data. This may lead to over-fitting! Performing cross-validation is recommended.');
    end
    uniqueChildren = fastsortedunique(sort(cat(1, allChildren{:})));
    
    for valItr = 1:validationFolds
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
        remainingChildren = fastsortedunique(sort(cat(1, allChildren{validSubIdx}))); %#ok<PFBNS>
        
        % If the maximum possible coverage has dropped, we should signal
        % that the matching threshold should be higher here.
        if (numel(remainingChildren) / numel(uniqueChildren)) < stoppingCoverage
           display(['[SUBDUE - Warning] We are losing our coverage on validation data: %' ...
               num2str(floor(100 * (numel(remainingChildren) / numel(uniqueChildren)))) ...
               ' maximal coverage. You may need to increase maximal threshold bound.']);
        end

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
       isSolutionOptimal = false;
       
       % For discriminative search, we need more than min, max and mid
       % values. Let's keep everything. 
       evaluatedThrs = [];
       evaluatedAccs = [];
       if supervisionFlag
           thrStack = [(minThr + midThr)/2, midThr, (midThr + maxThr)/2, maxThr];
           midThr = minThr;
       else
           thrStack = [];
       end

       display(['[SUBDUE] Starting threshold search with ' num2str(midThr) '. We are limited to ' num2str(numberOfFinalSubs) ' compositions.']);
       display(['[SUBDUE] Initially, we have ' num2str(numberOfBestSubs) ' subs to consider.']);
       %% Searching for optimal threshold using a binary search mechanism.
       while (currentDepth <= maxDepth) && ~isSolutionOptimal
           isCoverageOptimal = false;
           % If this is the final leg of search, we allow part selection to
           % have more subs than the max number.
           if currentDepth == maxDepth
              moreSubsAllowed = 1;
           end
           
           %% Optimality logic is implemented here. 
           % We're searching for an optimal threshold.
           % First, we obtain non-overlapping, minimal set of subs 
           % (common in both supervised and unsupervised learning)..
           [validSubs, overallCoverage, overallMatchScore] = getReconstructiveParts(bestSubs, ...
                numberOfFinalSubs, moreSubsAllowed, smartSubElimination, midThr, ...
                stoppingCoverage, remainingChildren, nodeDistanceMatrix, ...
                edgeDistanceMatrix, singlePrecision);
            
           % In order to measure the discrimination capabilities of our
           % system, we estimate the discrimination performance of the
           % selected subs. 
           accuracy = calculateCategorizationAccuracy(bestSubs(validSubs), categoryArrIdx, imageIdx, validationIdx, valItr, midThr, singlePrecision);
%          accuracy = calculateCategorizationAccuracy(bestSubs, categoryArrIdx, imageIdx, validationIdx, valItr, midThr, singlePrecision);
           
           %% We determine optimality of the solution based on the criterion.
           % If supervisionFlag is true, our criterion is classification
           % accuracy. Otherwise, we focus on coverage on the data.
           if supervisionFlag
               if accuracy > 1 - singlePrecision
                  isSolutionOptimal = true; 
               end
           else
               % If the coverage is good, mark the corresponding flag.
               if overallCoverage >= stoppingCoverage
                  isCoverageOptimal = true; 
               end
               
               % This is the multiplication to find our metric.
               unsupMetric = overallCoverage * overallMatchScore;

               % If the score is ideal, we mark solution flag as true.
               if unsupMetric >= stoppingCoverage
                   isSolutionOptimal = true;
               end
           end

           % Update final part list.
           optimalThreshold = midThr;

           %% Depending on the feedback, we're continuing or stopping the search.
            % Depending on which side we are on the ideal hidden node count, 
            % we make a binary search on the "good" threshold.
            currentDepth = currentDepth + 1;
            if (currentDepth <= maxDepth)
                if ~supervisionFlag
                    % We are performing unsupervised learning. Our
                    % threshold search is based on reconstruction error,
                    % not discrimination error.
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
                else
                    % We have switched to supervised learning already. The
                    % threshold search is driven by classification error.
                    evaluatedThrs = [evaluatedThrs; midThr]; %#ok<AGROW>
                    evaluatedAccs = [evaluatedAccs; accuracy]; %#ok<AGROW>
                    if ~isempty(thrStack)
                        midThr = thrStack(1);
                        thrStack = thrStack(2:end);
                    else
                        % Get maximum accuracy value for all evaluated thresholds (apart
                        % from both extremes), and select a new threshold
                        % which will reduce the uncertainty most.
                        [~, maxAccThrIdx] = max(evaluatedAccs);
                        differences = evaluatedThrs - evaluatedThrs(maxAccThrIdx);
                        posDifferences = differences;
                        negDifferences = differences;
                        posDifferences(posDifferences <= 0) = NaN;
                        negDifferences(negDifferences >= 0) = NaN;
                        
                        % Determine if we are further away from the
                        % next smaller threshold, or next greater
                        % threshold. Whicever is further, we take the
                        % mid-point between the local peak and it, and try
                        % it as our new threshold.
                        if isnan(min(posDifferences)) || min(posDifferences) < abs(max(negDifferences))
                           [~, minThrIdx] = max(negDifferences);
                           minThrVal = evaluatedThrs(minThrIdx);
                           midThr = (evaluatedThrs(maxAccThrIdx) + minThrVal) / 2;
                        else
                           [~, maxThrIdx] = min(posDifferences);
                           maxThrVal = evaluatedThrs(maxThrIdx);
                           midThr = (evaluatedThrs(maxAccThrIdx) + maxThrVal) / 2;
                        end
                    end
                end
            end
       end
       
       % If supervision flag is on, we simply select the best threshold
       % among all tried.
       if supervisionFlag
          [optimalAccuracy, maxThrIdx] = max(evaluatedAccs);
          if nnz(evaluatedAccs == optimalAccuracy) > 1
              optimalThreshold = min(evaluatedThrs(evaluatedAccs == optimalAccuracy));
          else
              optimalThreshold = evaluatedThrs(maxThrIdx);
          end
          crossValAccuracy(valItr) = optimalAccuracy;
       else
          crossValAccuracy(valItr) = accuracy;
       end
       crossValThresholds(valItr) = optimalThreshold;
       subEliminationFlags(valItr) = smartSubElimination;
    end
    
    %% We have found the optimal threshold. 
    % Now, we obtain the subs one final time using the new threshold, and exit.
    bestSubs = orgBestSubs;
    if validationFolds > 1
        display(['[SUBDUE] Thresholds learned from cross-validation folds: ' ...
            mat2str(crossValThresholds) ', with mean: ' num2str(median(crossValThresholds)) '.']);
    else
        display(['[SUBDUE] Matching threshold is determined as ' num2str(median(crossValThresholds))]); 
    end
    optimalThreshold = median(crossValThresholds);
    optimalAccuracy = mean(crossValAccuracy);
    smartSubElimination = mean(subEliminationFlags) > 0;
    uniqueChildren = fastsortedunique(sort(cat(1, allChildren{:})));
    display(['[SUBDUE] Mean cross-validation accuracy on the data: %' num2str(100 * optimalAccuracy) '.']); 
    
    % Get a new set of subs.
   [finalSubList, ~, ~] = getReconstructiveParts(bestSubs, ...
    numberOfFinalSubs, 1, smartSubElimination, optimalThreshold, ...
    stoppingCoverage, uniqueChildren, nodeDistanceMatrix, ...
    edgeDistanceMatrix, singlePrecision);

   % Update instance information.
   bestSubs = bestSubs(finalSubList);
   % Update bestSubs instances by taking the new threshold into account.
   parfor bestSubItr = 1:numel(bestSubs)
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