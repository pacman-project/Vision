%> Name: getDiscriminativeParts
%>
%> Description: This function is used to select a number of discriminative
%> parts from bestSubs. 
%>
%> @param bestSubs The list of selected subs so far. We will further
%> select a subset of this set, based on dicrimination.
%>
%> @retval validSubs Array that includes linear ids of selected subs in
%> bestSubs.
%> @retval overallCoverage The coverage of the data using the selected subs.
%> @retval overallMatchScore The average matching cost of parts in the
%> selected set, and the instances in the data.
%>
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 12.05.2015
function [validSubs, overallCoverage, overallMatchCost] = getDiscriminativeParts(bestSubs, numberOfFinalSubs, valItr, categoryArrIdx, ...
    imageIdx, validationIdx, smartSubElimination, midThr, optimalCount, uniqueChildren, nodeDistanceMatrix, ...
    edgeDistanceMatrix, singlePrecision)

   numberOfBestSubs = numel(bestSubs);
   remainingCategories = unique(categoryArrIdx);
   numberOfRemainingCategories = numel(remainingCategories);
   minProb = 1/numberOfRemainingCategories;
   prevGraphNodeCount = numel(uniqueChildren);
   
   % Calculate which leaf nodes are covered.
   subMatchCosts = zeros(numberOfBestSubs,1);
   subNodes = cell(numberOfBestSubs,1);
   parfor bestSubItr = 1:numberOfBestSubs
       adaptiveThreshold = ((midThr * (size(bestSubs(bestSubItr).edges,1) * 2 + 1)) + singlePrecision);
       validInstanceIdx = bestSubs(bestSubItr).instanceMatchCosts < adaptiveThreshold & bestSubs(bestSubItr).instanceValidationIdx ~= valItr;
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
       subMatchCosts(bestSubItr) = mean(matchCosts) / (size(bestSubs(bestSubItr).edges,1) * 2 + 1);
   end
   
   % Calculate which leaf nodes are covered.;
   subDiscriminativeScores = zeros(numberOfBestSubs,1);
   parfor bestSubItr = 1:numberOfBestSubs
       adaptiveThreshold = ((midThr * (size(bestSubs(bestSubItr).edges,1) * 2 + 1)) + singlePrecision);
       validInstanceIdx = bestSubs(bestSubItr).instanceMatchCosts < adaptiveThreshold & ...
           bestSubs(bestSubItr).instanceValidationIdx ~= valItr;
       instanceCategories = bestSubs(bestSubItr).instanceCategories(validInstanceIdx, :);
       
       % Get posterior probabilities of classes given the part, and obtain
       % the maximum probability. Score of the sub depends on it.
       posteriorProbs = hist(double(instanceCategories), double(remainingCategories));
       posteriorProbs = posteriorProbs / sum(posteriorProbs);
       maxProb = max(posteriorProbs);
       
       % Find sub score.
       subScore = (maxProb - minProb) * numel(instanceCategories);
 %      subScore = maxProb - minProb;
       subDiscriminativeScores(bestSubItr) = subScore;
   end

   % Eliminate subs that match to better subs. We'll have subs that are
   % far away from each other (in terms of pairwise distances), and
   % thus we will have less subs to evaluate.]
   if smartSubElimination
       validSubIdx = getDisjointSubs(bestSubs, ...
           nodeDistanceMatrix, edgeDistanceMatrix, singlePrecision, midThr);
       display(['[SUBDUE] [Val ' num2str(valItr) ']Considering only disjoint examples, we are down to ' num2str(nnz(validSubIdx)) ' subs to consider.']);
   else
       validSubIdx = ones(numel(bestSubs),1) > 0;
   end
   validSubIdx = find(validSubIdx);
   
   %% Get top N discriminative subs.
   [subScores, orderedSubIdx] = sort(subDiscriminativeScores(validSubIdx), 'descend');
   subScores = subScores/max(subScores);
   
   % Greedily select best subs.
   bestAcc = 0;
   valAccuracyArr = zeros(numel(orderedSubIdx),1);
   stepSize = 20;
   if optimalCount == -1 && numel(orderedSubIdx) >= 100
       subCheckPoints = stepSize:stepSize:min(numel(orderedSubIdx), numberOfFinalSubs);
       for bestSubItr = stepSize:stepSize:min(numel(orderedSubIdx), numberOfFinalSubs)
           [accuracy, ~] = calculateCategorizationAccuracy(bestSubs(validSubIdx(orderedSubIdx(1:bestSubItr))), ...
               categoryArrIdx, imageIdx, validationIdx, valItr, midThr, singlePrecision, 1, false);
           valAccuracyArr(bestSubItr) = accuracy;
           if bestAcc < accuracy
               bestAcc = accuracy;
           end
       end
       
       % We select a cutting point.
       smoothingMask = makeGauss1D(1);
       smoothingMask = smoothingMask(2:(end-1));
       valAccuracyArrSmooth = convolve1D(valAccuracyArr(stepSize:stepSize:min(numel(orderedSubIdx), numberOfFinalSubs)), smoothingMask);
       [~, maxLoc] = max(valAccuracyArrSmooth);
       maxLoc = maxLoc + 2;
       addedSubs = 1:subCheckPoints(maxLoc);
       % Finally, select the subs.
   elseif optimalCount > -1
       addedSubs = 1:min(optimalCount, numel(orderedSubIdx));
   else
       addedSubs = 1:numel(orderedSubIdx);
   end
   
    % Get train and validation accuracy for final evaluation.
   [trainAccuracy, ~] = calculateCategorizationAccuracy(bestSubs(validSubIdx(orderedSubIdx(addedSubs))), ...
       categoryArrIdx, imageIdx, validationIdx, valItr, midThr, singlePrecision, 0, true);
   [trueAccuracy, ~] = calculateCategorizationAccuracy(bestSubs(validSubIdx(orderedSubIdx(addedSubs))), ...
       categoryArrIdx, imageIdx, validationIdx, valItr, midThr, singlePrecision, 1, true);
   display(['[SUBDUE] [Val ' num2str(valItr) '] Val accuracy after discriminative part selection : %' num2str(100 * trueAccuracy) ...
        ', with training set accuracy : %' num2str(100 * trainAccuracy), ...
        ', having ' num2str(numel(addedSubs)) ' subs.']);    
   validSubs = sort(validSubIdx(orderedSubIdx(addedSubs)));
   
   if isempty(validSubs)
       validSubs = 1:numel(subScores);
   end
   
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
   overallMatchCost = mean(subMatchCosts(validSubs));
   
   % Printing.
   display(['[SUBDUE] [Val ' num2str(valItr) '] We have selected  ' num2str(numel(validSubs)) ...
        ' out of ' num2str(numberOfBestSubs) ' subs.. Coverage: ' num2str(overallCoverage) ', average normalized match cost:' num2str(overallMatchCost) '.']);

end