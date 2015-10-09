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
            realEdgeLabels, nodePositions, edgeCoords, allEdges, allEdgeProbs, numberOfFinalSubs, stoppingCoverage, uniqueChildren)

   minNodeProbability = 0.00001;
   numberOfBestSubs = numel(bestSubs);
   prevGraphNodeCount = numel(uniqueChildren);
   maxChildId = max(uniqueChildren);
   prevGraphNodeLogProbs = zeros(1, maxChildId, 'single');
   minLogProb = single(log(minNodeProbability));
   allEdgesArr = {allEdges.adjInfo};
   nonzeroIdx = cellfun(@(x) ~isempty(x), allEdgesArr);
   allEdgesPeripheralNodes = cell(size(allEdgesArr));
   allEdgesPeripheralNodes(nonzeroIdx) = cellfun(@(x) x(:,2), allEdgesArr(nonzeroIdx), 'UniformOutput', false);
   
   %% Initially, we learn distributions of data points given each sub's center.
   RFSize = sqrt(size(edgeCoords,1));
   largeRFSize = RFSize * 2 - 1;
   largeRFSizes = [largeRFSize, largeRFSize];
   halfSize = RFSize;
   
   % Allocate space for node/edge distributions, both discrete.
   % We start by calculating number of sub, child pairs.
   maxLabel = max(realNodeLabels);
   numberOfSubChildPairs = 0;
   for subItr = 1:numberOfBestSubs
       numberOfSubChildPairs = numberOfSubChildPairs + size(bestSubs(subItr).edges,1) + 1;
   end
   
   % TODO: The following two entries can be made sparse if we're short on memory.
   labelProbArr = zeros(numberOfSubChildPairs, maxLabel);
   posProbArr = zeros(numberOfSubChildPairs, largeRFSize, largeRFSize);
   subInstancePositions = cell(numberOfBestSubs,1);
   pairCounter = 1;
   for subItr = 1:numberOfBestSubs
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
       subInstancePositions{subItr} = instancePositions;
       
       % Finally, collect statistics for every child.
       for childItr = 1:numberOfChildren
          % Learn node label distribution
          nodeLabels = double(realNodeLabels(instanceChildren(:,childItr)));
          entries = unique(nodeLabels);
          if numel(entries) > 1
              [nodeProbs, ~] = hist(nodeLabels, entries);
              nodeProbs = nodeProbs / sum(nodeProbs);
          else
              nodeProbs = 1;
          end
          labelProbArr(pairCounter, entries) = nodeProbs;
          
           % Learn position distributions
           % TODO: Make the distributions more continuous, as in gaussians.
           % Right now, they're entirely discrete.
           relativePositions = (nodePositions(instanceChildren(:,childItr),:) - instancePositions) + halfSize;
           relativePositionIdx = double(sub2ind(largeRFSizes, ...
               relativePositions(:,1), relativePositions(:,2)));
           uniquePosIdx = double(unique(relativePositionIdx));
           if numel(uniquePosIdx) > 1
              [posProbs, ~] = hist(relativePositionIdx, uniquePosIdx);
           else
               posProbs = 1;
           end
           posProbs = posProbs / sum(posProbs);
           probSlice = squeeze(posProbArr(pairCounter, :,:));
           probSlice(uniquePosIdx) = posProbs;
           posProbArr(pairCounter, :,:) = probSlice;
           
           % Increase counter.
           pairCounter = pairCounter + 1;
       end
   end
   
   % Check for a potential error condition.
   if size(posProbArr,2) ~= largeRFSize || size(posProbArr,3) ~= largeRFSize
      error('Problem in getReconstructiveParts: Relative coordinations are wrong!'); 
   end
   
   %% Start the main loop. We select subs one by one.
   % This is a greedy selection algorithm. Each new sub's evaluation
   % depends on the previous set, i.e. it shows how much novelty it can bring to the table.
   curLabelLogProbs = prevGraphNodeLogProbs;
   curLabelProbPairs = ones(maxChildId, maxLabel);
   curPosLogProbs = curLabelLogProbs;
   curPosProbPairs = sparse(zeros(maxChildId, 1));
   coveredNodes = zeros(1, maxChildId) > 0;
   selectedSubIdx = zeros(numberOfBestSubs,1) > 0;
   selectedSubs = [];
   addedSubOffset = 1;   
   pairCounter = 1;
   reconstructionArr = zeros(1, maxChildId) > 0;
   while addedSubOffset <= numberOfFinalSubs
       % Select next best sub.
       valueArr = -inf(numberOfBestSubs,1);
       curValue = -inf;
       for bestSubItr = 1:numberOfBestSubs
            % If this sub has already been selected, we move on.
            if selectedSubIdx(bestSubItr)
               continue;
            end
            
            instanceChildren = bestSubs(bestSubItr).instanceChildren;
            numberOfInstances = size(instanceChildren,1);
            allAssgnProbs = [];
            for childItr = 1:size(instanceChildren,2)
                allAssgnProbs = [allAssgnProbs; ...
                    repmat(labelProbArr(pairCounter,:), numberOfInstances, 1)]; %#ok<AGROW>
            
                % Increase counter.
                pairCounter = pairCounter + 1;
            end
            instanceChildrenRep = instanceChildren(:);
            
            % Obtain unique numbers.
            [uniqueChildren, IA, IC] = unique(instanceChildrenRep, 'stable');
            for childItr = 1:numel(IA)
               curLabelProbPairs(uniqueChildren(childItr),:) = prod(cat(1, curLabelProbPairs(uniqueChildren(childItr),:), ...
                   allAssgnProbs(IC == childItr, :)), 1);
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
   dataLikelihood = abs(mean(curLogProbs)/ minLogProb);
   
   % Record preserved subs.
   validSubs = find(selectedSubIdx);
   
   % Printing.
   display(['[SUBDUE] We have selected  ' num2str(numel(validSubs)) ...
        ' out of ' num2str(numberOfBestSubs) ' subs.. Coverage: ' num2str(overallCoverage) ', average normalized match cost:' num2str(overallMatchCost) '.']);

end