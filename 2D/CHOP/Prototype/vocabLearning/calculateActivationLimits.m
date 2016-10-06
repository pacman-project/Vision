%> Name: calculateActivationLimits
%>
%> Description: For every part, this function looks at each instance, and
%> calculates the minimal activation in the training set. The algorithm for
%> this method is to cover the whole dataset with maximal activation, and
%> then eliminate low-activation realizations that are deemed unnecessary.
%> 
%> @param vocabLevel Vocabulary layer.
%> @param graphLevel Realizations.
%> @param numberOfLeafNodes Number of leaf nodes to cover.
%>
%> @retval vocabLevel Parts with updated thresholds.
%> @retval graphLevel Part realizations (remaining).
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 05.10.2016
function [ vocabLevel, graphLevel ] = calculateActivationLimits( vocabLevel, graphLevel, numberOfLeafNodes )
     activations = cat(1, graphLevel.activation);
     labelIds = cat(1, graphLevel.realLabelId);
     leafNodes = {graphLevel.leafNodes};
     leafNodeActivationArr = -inf(numberOfLeafNodes,1);
     leafNodeLabelArr = zeros(numberOfLeafNodes,1);
     numberOfInstances = numel(graphLevel);
     
     %% Using a greedy approach, we lower the thresholds.
     for instanceItr = 1:numberOfInstances
          tempNodes  = leafNodes{instanceItr};
          tempLabels = leafNodeLabelArr(tempNodes);
          tempValsOld = leafNodeActivationArr(tempNodes);
          
          tempVals = max(tempValsOld, activations(instanceItr));
          
          if ~isequal(tempValsOld,tempVals)
               tempLabels(tempVals > tempValsOld) = labelIds(instanceItr);
               leafNodeLabelArr(tempNodes) = tempLabels;
               leafNodeActivationArr(tempNodes) = tempVals;
          end
     end
     
     %% At this step, we fit a gaussian to the data, and try to contain %99 of the samples by putting a new limit.
     for vocabNodeItr = 1:numel(vocabLevel)
          % Obtain data.
          dataToFit = activations(labelIds == vocabNodeItr);
          
          % Replicate the data by taking a mirror image (Make it a bell,
          % rather than half bell).
          [bins, binCenters] = hist(dataToFit, 10);
          [~, maxIdx] = max(bins);
          centerPoint = binCenters(maxIdx);
          dataToFit = cat(1,-dataToFit + 2 * centerPoint, dataToFit);
          
          % Fit gaussian, take %99.7 interval.
          mu = mean(dataToFit);
          var = std(dataToFit);
          vocabLevel(vocabNodeItr).minActivationLog = max(vocabLevel(vocabNodeItr).minActivationLog, mu - var * 2.5);
     end
     
     %% Go through the list and record activation threshold.
     validNodeArr = ones(numberOfInstances,1) > 0;
     for vocabNodeItr = 1:numel(vocabLevel)
          minThr = max(vocabLevel(vocabNodeItr).minActivationLog, min(leafNodeActivationArr(leafNodeLabelArr == vocabNodeItr)));
          if isempty(minThr)
               minThr = 0;
          end
          vocabLevel(vocabNodeItr).minActivationLog = single(minThr);
          validNodeArr(labelIds == vocabNodeItr) = validNodeArr(labelIds == vocabNodeItr) & activations(labelIds == vocabNodeItr) >= single(minThr);
     end
     
     %% Remove activations that do not pass the threshold values.
     graphLevel = graphLevel(validNodeArr);
end