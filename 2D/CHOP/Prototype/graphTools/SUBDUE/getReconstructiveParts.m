%> Name: getReconstructiveParts
%>
%> Description: Given a set of parts in bestSubs, this function greedily
%> selects a set of parts that minimize the likelihood of the data. The data
%> is grouped into overlapping receptive fields, and the reduction in the
%> cost is associated with increasing likelihood of the underlying data. Two
%> factors are contributing towards the data likelihood description, namely
%> node label and position prediction. 
%> 
%> @param bestSubs Initial set of substructures.
%> @param realNodeLabels The real labels of underlying data. 
%> @param realEdgeLabels The real labels of the edges that encode spatial
%> distributions in the bottom level. 
%> @param allEdges All edges encoded in the first level, with each cell
%> corresponding to a separate node's edges. 
%> @param allEdgeProbs Probabilities associated with edges.
%> @param numberOfFinalSubs Selection will stop if the number of selected
%> subs exceeds numberOfFinalSubs. 
%> @param stoppingCoverage The minimum coverage that is required to stop
%> selection.
%> @param uniqueChildren The ids of the nodes(data) to be covered.
%>
%> @retval validSubs Ids of final subs.
%> @retval overallCoverage The coverage on the data.
%> @retval dataLikelihood Data likelihood given the selected parts.
%>
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 08.10.2015
function [validSubs, overallCoverage, dataLikelihood] = getReconstructiveParts(bestSubs, realNodeLabels, ...
            nodePositions, edgeCoords, numberOfFinalSubs, uniqueChildren)

   minNodeProbability = 0.00001;
   coverageStoppingVal = 0.999;
   likelihoodStoppingVal = 0.99;
   coveragePartSelectionMethod = 'Greedy'; % 'Genetic' or 'Greedy'
   % Override part selection method if we have too many parts to select
   % from.
   if numberOfBestSubs > 2000
        coveragePartSelectionMethod = 'Greedy';
   end
   numberOfBestSubs = numel(bestSubs);
   prevGraphNodeCount = numel(uniqueChildren);
   maxChildId = max(uniqueChildren);
   prevGraphNodeLogProbs = zeros(1, maxChildId, 'single');
   minLogProb = single(log(minNodeProbability));
   
   %% Initially, we learn distributions of data points given each sub's center.
   RFSize = sqrt(size(edgeCoords,1));
   largeRFSize = RFSize * 2 - 1;
   largeRFSizes = [largeRFSize, largeRFSize];
   halfSize = RFSize;
   
   % Allocate space for node/edge distributions, both discrete.
   % We start by calculating number of sub, child pairs.
   numberOfSubChildPairs = 0;
   for subItr = 1:numberOfBestSubs
       numberOfSubChildPairs = numberOfSubChildPairs + size(bestSubs(subItr).edges,1) + 1;
   end
   
   % Allocate space for log likelihood results.
   posProbArr = zeros(numberOfSubChildPairs, largeRFSize, largeRFSize);
   subLabelProbs = cell(numberOfBestSubs, 1);
   subPosProbs = cell(numberOfBestSubs, 1);
   subCoveredNodes = cell(numberOfBestSubs,1);
   
   % Go over all possible part-subpart pairs, and calculate probabilities.
   parfor subItr = 1:numberOfBestSubs
       instanceChildren = bestSubs(subItr).instanceChildren;
       numberOfInstances = size(instanceChildren,1);
       numberOfChildren = size(instanceChildren,2);
       instanceMappings = bestSubs(subItr).instanceMappings;
       instancePositions = zeros(numberOfInstances, 2, 'int32');
       
       % If the mappings are not trivial (1:numberOfChildren at every row),
       % we re-order children to reflect the mappings.
       if ~issorted(instanceMappings, 'rows');
           for instanceItr = 1:numberOfInstances
              instanceChildren(instanceItr,:) = instanceChildren(instanceItr, ...
                  instanceMappings(instanceItr,:));
           end
       end
       bestSubs(subItr).instanceChildren = instanceChildren;
       
       % Find center positions for each instance.
       for instanceItr = 1:numberOfInstances
          instancePositions(instanceItr,:) = int32(round(sum(nodePositions(instanceChildren(instanceItr,:), :),1) ...
                               / numberOfChildren)); 
       end
       labelLogProbs = prevGraphNodeLogProbs;
       posLogProbs = labelLogProbs;
       
       % Finally, collect statistics for every child.
       for childItr = 1:numberOfChildren
          % Learn node label distribution
          nodeLabels = double(realNodeLabels(instanceChildren(:,childItr))); %#ok<*PFBNS>
          entries = unique(nodeLabels);
          inverseArr = zeros(1, max(entries));
          inverseArr(entries) = 1:numel(entries);
          if numel(entries) > 1
              [nodeProbs, ~] = hist(nodeLabels, entries);
              nodeProbs = nodeProbs / sum(nodeProbs);
          else
              nodeProbs = 1;
          end
          labelLogProbs(instanceChildren(:,childItr)) = max(labelLogProbs(instanceChildren(:,childItr)), log(nodeProbs(inverseArr(nodeLabels))) - minLogProb);
          
           % Learn position distributions
           % TODO: Make the distributions more continuous, as in gaussians.
           % Right now, they're entirely discrete.
           relativePositions = (nodePositions(instanceChildren(:,childItr),:) - instancePositions) + halfSize;
           relativePositionIdx = double(sub2ind(largeRFSizes, ...
               relativePositions(:,1), relativePositions(:,2)));
           uniquePosIdx = double(unique(relativePositionIdx));
           inverseArr = zeros(1, max(uniquePosIdx));
           inverseArr(uniquePosIdx) = 1:numel(uniquePosIdx);
           if numel(uniquePosIdx) > 1
              [posProbs, ~] = hist(relativePositionIdx, uniquePosIdx);
           else
               posProbs = 1;
           end
           posProbs = posProbs / sum(posProbs);
           
           % Assign position probabilities to samples, and save the info.
           posLogProbs(instanceChildren(:,childItr)) = max(posLogProbs(instanceChildren(:,childItr)), log(posProbs(inverseArr(relativePositionIdx))) - minLogProb);
       end
        % Save child probabilities.
        allNodes = unique(instanceChildren);
        allNodes = (allNodes(:))';
        subCoveredNodes{subItr} = allNodes;
        subLabelProbs{subItr} = labelLogProbs(allNodes);
        subPosProbs{subItr} = posLogProbs(allNodes);
   end
   
   % Check for a potential error condition.
   if size(posProbArr,2) ~= largeRFSize || size(posProbArr,3) ~= largeRFSize
      error('Problem in getReconstructiveParts: Relative coordinations are wrong!'); 
   end
   
   % Create a genetic algorithm instance to solve the part selection
   % problem.
   f = @(x) paramfunc(x, subLabelProbs, subPosProbs, subCoveredNodes);
   g = @(x) paramfunc2(x, subLabelProbs, subPosProbs, subCoveredNodes);
   h = @(x) paramfunc3(x, subLabelProbs, subPosProbs, subCoveredNodes);
     
     %% Finally, we implement an algorithm for part selection. 
     % This step is done to perform an initial pass to reduce the number of
     % parameters (subs to be selected) significiantly. Then, we perform
     % another pass using a data likelihood measure. This step implements a
     % coverage-based part selection mechanism.
     
     if strcmp(coveragePartSelectionMethod, 'Genetic')
          A = ones(1, numberOfBestSubs);
          b = numberOfFinalSubs;
          LB = zeros(1,numberOfBestSubs);
          UB = ones(1,numberOfBestSubs);
          IntCon = 1:numberOfBestSubs;
          %   IntCon = [];
          selectedSubIdx = ones(1,numberOfBestSubs) > 0;
          stoppingFVal = coverageStoppingVal * f(selectedSubIdx);
          options = gaoptimset('Display', 'diagnose', 'UseParallel', 'Always', 'FitnessLimit', stoppingFVal, 'Generations', 1000, ...
               'CreationFcn', @gacreationlinearfeasible, 'CrossoverFcn', @crossoverintermediate, 'HybridFcn', @fminsearch, ...
               'MutationFcn', @mutationadaptfeasible, 'PopulationSize', 1000, 'TolFun', 1e-8);
          [selectedSubIdx, fval, exitFlag] = ga(h, numberOfBestSubs, A, b, [], [], LB, UB, [], IntCon, options);
          display(['Exit flag for genetic algorithm: ' num2str(exitFlag) ', with fval: ' num2str(fval) '.']);
     else
          subCounter = 0; 
          selectedSubIdx = zeros(1,numberOfBestSubs) > 0;
          curFVal = 0;
          valueArr = inf(1,numberOfBestSubs);
          while subCounter < numberOfFinalSubs
               maxLocVal = 0;
               maxSubIdx = [];

               for subItr = 1:numberOfBestSubs
                    if selectedSubIdx(subItr) == 1 || maxLocVal > valueArr(subItr)
                         continue;
                    end
                    tempSubIdx = selectedSubIdx;
                    tempSubIdx(subItr) = 1;
                    diffVal = -h(tempSubIdx);

                    % Save diffVal.
                    if diffVal < valueArr(subItr)
                         valueArr(subItr) = diffVal;
                    end

                    % Save value if this part has maximum value.
                    if diffVal > 0 && diffVal > maxLocVal
                         maxLocVal = diffVal;
                         maxSubIdx = tempSubIdx;
                    end
               end
               valueArr(maxSubIdx) = 0;
               if isempty(maxSubIdx)
                    break;
               end
               curFVal = h(maxSubIdx);
               selectedSubIdx = maxSubIdx;
               subCounter = subCounter + 1;
          end
     end
   fval = curFVal;
   validSubs = find(selectedSubIdx >= 0.5);
   subCoveredNodes = subCoveredNodes(validSubs);
   subLabelProbs = subLabelProbs(validSubs);
   subPosProbs = subPosProbs(validSubs);
   
     %% Now, it's time for a data likelihood-driven part selection.
   %    A = ones(1, numberOfBestSubs);
%    b = numberOfFinalSubs;
%    LB = zeros(1,numberOfBestSubs);
%    UB = ones(1,numberOfBestSubs);
% %   IntCon = 1:numberOfBestSubs;
%    IntCon = [];
%    selectedSubIdx = ones(1,numberOfBestSubs) > 0;
%    stoppingFVal = coverageStoppingVal * f(selectedSubIdx);
%    options = gaoptimset('Display', 'diagnose', 'UseParallel', 'Always', 'FitnessLimit', stoppingFVal, 'Generations', 1000, ...
%    'CreationFcn', @gacreationlinearfeasible, 'CrossoverFcn', @crossoverintermediate, 'HybridFcn', @fminsearch, ...
%    'MutationFcn', @mutationadaptfeasible, 'PopulationSize', 1000, 'TolFun', 1e-8);
%    [selectedSubIdx, fval, exitFlag] = ga(h, numberOfBestSubs, A, b, [], [], LB, UB, [], IntCon, options);
   
   % Calculate statistics.
   allCoveredNodes = cat(2, subCoveredNodes{:});
   allCoveredNodes = unique(allCoveredNodes);
   overallCoverage = numel(allCoveredNodes) / prevGraphNodeCount;
%   dataLikelihood = 1 - (fval / prevGraphNodeCount - minLogProb * 2) / (-minLogProb * 2);
   dataLikelihood = fval/prevGraphNodeCount;
   
   % Printing.
   display(['[SUBDUE] We have selected  ' num2str(numel(validSubs)) ...
        ' out of ' num2str(numberOfBestSubs) ' subs.. Coverage: ' num2str(overallCoverage) ', average data point likelihood:' num2str(dataLikelihood) '.']);

end

% Fitness function 1 for the genetic algorithm.
function y = paramfunc(x, subLabelProbs, subPosProbs, subCoveredNodes)   
     x = x>= 0.5;
     allSubLabelProbs = cat(2, subLabelProbs{x});
     allSubPosProbs = cat(2, subPosProbs{x});
     allCoveredNodes = cat(2, subCoveredNodes{x});
     
     %% Get max label/pos log probs for every child.
     % First, we obtain max label log probs.
     [labelVals, labelSortIdx] = sort(allSubLabelProbs, 'descend');
     sortedCoveredNodes = allCoveredNodes(labelSortIdx);
     [~, IA, ~] = unique(sortedCoveredNodes, 'first');
     r1 = labelVals(IA);
     
     %Then, we obtain max location log probs.
     [posVals, posSortIdx] = sort(allSubPosProbs, 'descend');
     sortedCoveredNodes = allCoveredNodes(posSortIdx);
     [~, IA, ~] = unique(sortedCoveredNodes, 'first');
     r2 = posVals(IA);
     
     y = -round(sum(double(r1)) + sum(double(r2)));
end

% Fitness function 2 for the genetic algorithm.
function y = paramfunc2(x, subLabelProbs, subPosProbs, subCoveredNodes)   
     x = x>=0.5;
     allSubLabelProbs = cat(2, subLabelProbs{x});
     allSubPosProbs = cat(2, subPosProbs{x});
     allCoveredNodes = cat(2, subCoveredNodes{x});
     
     % The following objective function implementation does not
     % differentiate label, position probability contributions from each
     % prediction.
    logProbs = allSubLabelProbs + allSubPosProbs;
    
     %% Get max log probs for every child.
     [vals, idx] = sort(logProbs, 'descend');
     sortedCoveredNodes = allCoveredNodes(idx);
     [~, IA, ~] = unique(sortedCoveredNodes, 'stable');
     vals = vals(IA);
     y = round(-sum(double(vals)));
end

% Fitness function 3 for the genetic algorithm.
function y = paramfunc3(x, ~, ~, subCoveredNodes)   
     x = x>=0.5;
     allCoveredNodes = cat(2, subCoveredNodes{x});
     numberOfNodes = numel(allCoveredNodes);
     allCoveredNodes = fastsortedunique(sort(allCoveredNodes));
     y = -2*numel(allCoveredNodes) + numberOfNodes;
end
