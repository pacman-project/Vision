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
function [ vocabLevel, graphLevel ] = calculateActivations(vocabLevel, vocabLevelDistributions, graphLevel, prevActivations, prevPosition, levelItr, options)
    rfSize = getRFSize(options, levelItr);
    halfRFSize = round(rfSize(1)/2);    
    maxPosProb = 0.1;

    %% Go through the list of vocabulary realizations and calculate activations.
    newActivations = zeros(numel(graphLevel), 1, 'single');
    graphLabels = [graphLevel.labelId];
    graphLevelChildren = {graphLevel.children};
    graphLevelPrecisePositions = cat(1, graphLevel.precisePosition);
    for vocabLevelItr = 1:numel(vocabLevel)
        %% First, we start by calculating position likelihood.
        idx = graphLabels == vocabLevelItr;
        instanceChildren = cat(1, graphLevelChildren{idx});
        instancePositions = graphLevelPrecisePositions(idx, :);
        instanceChildrenCombinedPos = zeros(size(instanceChildren,1), size(instanceChildren,2) * 2);
        
        % Put children positions in an array.
        for instanceItr = 1:size(instanceChildren,2)
             relevantChildren = instanceChildren(:, instanceItr);
             
             % If we're working with peripheral sub-parts, we calculate
             % position distributions as well.
             samples = prevPosition(relevantChildren, :) - instancePositions;

             % Save samples (positions and labels).
             instanceChildrenCombinedPos(:, ((instanceItr-1)*2+1):((instanceItr)*2)) = samples;
        end
        
        % Calculate probability.
        relevantSamples = instanceChildrenCombinedPos / halfRFSize;
        
        % Finally, we choose a pdf threshold for this sub, based on
        % the values of individual instances.
        obj = vocabLevelDistributions(vocabLevelItr).childrenPosDistributions;
        posActivations = pdf(obj, relevantSamples);

        % Convert to log probabilities
        activationDenom = max(pdf(obj, obj.mu)) * 1/maxPosProb;
        posActivations = posActivations / activationDenom;
        posActivations = single(log(posActivations));
        
        %% Second, we take the children's activations into account as well.
        prevChildrenActivations = prevActivations(instanceChildren);
        if ~isequal(size(prevChildrenActivations), size(instanceChildren))
            prevChildrenActivations = prevChildrenActivations';
        end
        avgActivations = mean([prevChildrenActivations, posActivations], 2);
        newActivations(idx) = avgActivations;
        
        % Save the minimum log value for further filtering.
        vocabLevel(vocabLevelItr).minActivationLog = min(avgActivations);
    end
    newActivations = num2cell(newActivations);
    [graphLevel.activation] = deal(newActivations{:});
end

