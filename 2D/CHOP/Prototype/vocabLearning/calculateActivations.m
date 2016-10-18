%> Name: calculateActivations
%>
%> Description: Calculates activation values for the data and given
%> realizations. For each realization, the level 1 parts which are not
%> explained are penalized with a small probability.
%> 
%> @param graphLevel List of realizations to calculate activations of. 
%> @param level1Coords Level 1 nodes, with reduced coordinates. 
%>
%> @retval extendedSubs Extended sub list.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 18.01.2016
function [ vocabLevel, vocabLevelDistributions, graphLevel ] = calculateActivations(vocabLevel, vocabLevelDistributions, graphLevel, prevActivations, prevPosition, levelItr, options)
    rfSize = getRFSize(options, levelItr);
    halfRFSize = round(rfSize(1)/2);
    realRFSize = halfRFSize * 2 + 1;
    minLogProb = single(-20);
    minProb = single(exp(-20));

    %% Go through the list of vocabulary realizations and calculate activations.
    newActivations = zeros(numel(graphLevel), 1, 'single');
    graphLabels = [graphLevel.labelId];
    graphLevelChildren = {graphLevel.children};
    graphLevelPrecisePositions = cat(1, graphLevel.precisePosition);
    newActivationArr = cell(numel(vocabLevel),1);
    idxArr = cell(numel(vocabLevel),1);
    
    for vocabLevelItr = 1:numel(vocabLevel)
        %% First, we start by calculating position likelihood.
        idx = graphLabels == vocabLevelItr;
        instanceChildren = cat(1, graphLevelChildren{idx});
        instancePositions = graphLevelPrecisePositions(idx, :);
        childrenProbs = vocabLevelDistributions(vocabLevelItr).childrenPosDistributionProbs;
        childrenPosActivations = zeros(size(instanceChildren), 'single');
        
        % Put children positions in an array.
        minPosProbs = zeros(size(instanceChildren,2),1);
        for instanceItr = 1:size(instanceChildren,2)
             relevantChildren = instanceChildren(:, instanceItr);
             relevantProbs = childrenProbs{instanceItr};
             minPosProbs(instanceItr) = log(min(min(relevantProbs(relevantProbs > 0))));
             
             % Calculate position probabilities.
             samples = (prevPosition(relevantChildren, :) - instancePositions) + halfRFSize + 1;
             validIdx = samples > 0 & samples < realRFSize + 1;
             validIdx = validIdx(:,1) & validIdx(:,2);
             pointProbs = ones(size(samples,1),1, 'single') * minProb;
             samplesIdx = samples(validIdx,1) + (samples(validIdx,2)-1)*realRFSize;
             pointProbs(validIdx) = full(relevantProbs(samplesIdx));
             childrenPosActivations(:, instanceItr) = log(pointProbs);
        end
        
        %% Second, we take the children's activations into account as well.
        prevChildrenActivations = prevActivations(instanceChildren);
        if ~isequal(size(prevChildrenActivations), size(instanceChildren))
            prevChildrenActivations = prevChildrenActivations';
        end
        if levelItr > 2
              avgActivations = mean([prevChildrenActivations, childrenPosActivations], 2);
        else
              avgActivations = mean(childrenPosActivations,2);
        end
        avgActivations(avgActivations < minLogProb) = minLogProb;
        
        % Calculate minimum possible probabilities.
        minPosLog = max(minLogProb, mean(minPosProbs));
        if size(childrenPosActivations,1) < 3
             minLog = mean([min(mean(prevChildrenActivations,2)), minPosLog]);
        else
             minLog = min(avgActivations);
        end
        
        vocabLevelDistributions(vocabLevelItr).minPosActivationLog = minPosLog;
        newActivationArr{vocabLevelItr} = avgActivations;
        idxArr{vocabLevelItr} = idx;
        
        % Save the minimum log value for further filtering.
        vocabLevel(vocabLevelItr).minActivationLog = minLog;
    end
    
    % Finally, assign activations back.
    for vocabLevelItr = 1:numel(vocabLevel)
         newActivations(idxArr{vocabLevelItr}) = newActivationArr{vocabLevelItr};
    end
    newActivations = num2cell(newActivations);
    [graphLevel.activation] = deal(newActivations{:});
end

