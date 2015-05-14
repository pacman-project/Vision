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
function [validSubs, overallCoverage, overallMatchScore] = getDiscriminativeParts(bestSubs, valItr, categoryArrIdx, ...
    numberOfFinalSubs, smartSubElimination, midThr, uniqueChildren, nodeDistanceMatrix, ...
    edgeDistanceMatrix, singlePrecision)

   numberOfBestSubs = numel(bestSubs);
   remainingCategories = unique(categoryArrIdx);
   numberOfRemainingCategories = numel(remainingCategories);
   minProb = 1/numberOfRemainingCategories;
   prevGraphNodeCount = numel(uniqueChildren);
   
   % Calculate which leaf nodes are covered.
   subNodes = cell(numberOfBestSubs,1);
   subMatchScores = cell(numberOfBestSubs,1);
   parfor bestSubItr = 1:numberOfBestSubs
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
   
   % Calculate which leaf nodes are covered.;
   subDiscriminativeScores = zeros(numberOfBestSubs,1);
   parfor bestSubItr = 1:numberOfBestSubs
       adaptiveThreshold = ((midThr * (size(bestSubs(bestSubItr).edges,1) * 2 + 1)) + singlePrecision);
       validInstanceIdx = bestSubs(bestSubItr).instanceMatchCosts < adaptiveThreshold;
       instanceCategories = bestSubs(bestSubItr).instanceCategories(validInstanceIdx, :);
       
       % Get posterior probabilities of classes given the part, and obtain
       % the maximum probability. Score of the sub depends on it.
       posteriorProbs = hist(double(instanceCategories), double(remainingCategories));
       posteriorProbs = posteriorProbs / sum(posteriorProbs);
       maxProb = max(posteriorProbs);
       
       % Find sub score.
       subScore = (maxProb - minProb) * numel(instanceCategories);
       subDiscriminativeScores(bestSubItr) = subScore;
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
   validSubIdx = find(validSubIdx);
   
   %% Get top N discriminative subs.
   [subScores, orderedSubIdx] = sort(subDiscriminativeScores(validSubIdx), 'descend');
   subScores = subScores/max(subScores);
   thr = graythresh(subScores);
   orderedSubIdx = orderedSubIdx(subScores >= thr);
%    if numel(orderedSubIdx)>numberOfFinalSubs
%        orderedSubIdx = orderedSubIdx(1:numberOfFinalSubs);
%    end
   validSubs = sort(validSubIdx(orderedSubIdx));
   
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
   display(['[SUBDUE] We have selected  ' num2str(numel(validSubs)) ...
        ' out of ' num2str(numberOfBestSubs) ' subs.. Coverage: ' num2str(overallCoverage) ', average match score:' num2str(overallMatchScore) '.']);

end