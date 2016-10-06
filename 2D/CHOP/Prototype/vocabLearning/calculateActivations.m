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
    halfPixel = 1/(2*halfRFSize);
    minLogProb = single(-20);
    minProb = single(exp(-20));

    %% Go through the list of vocabulary realizations and calculate activations.
    newActivations = zeros(numel(graphLevel), 1, 'single');
    graphLabels = [graphLevel.labelId];
    graphLevelChildren = {graphLevel.children};
    graphLevelPrecisePositions = cat(1, graphLevel.precisePosition);
    newActivationArr = cell(numel(vocabLevel),1);
    idxArr = cell(numel(vocabLevel),1);
    
    parfor vocabLevelItr = 1:numel(vocabLevel)
        %% First, we start by calculating position likelihood.
        idx = graphLabels == vocabLevelItr;
        instanceChildren = cat(1, graphLevelChildren{idx});
        instancePositions = graphLevelPrecisePositions(idx, :);
        relevantSamples = zeros(size(instanceChildren,1), size(instanceChildren,2) * 2);
        obj = vocabLevelDistributions(vocabLevelItr).childrenPosDistributions;
        mu = obj.mu;
        Sigma = obj.Sigma;
        childrenPosActivations = zeros(size(instanceChildren), 'single');
        
        % Put children positions in an array.
        for instanceItr = 1:size(instanceChildren,2)
             relevantChildren = instanceChildren(:, instanceItr);
             
             % If we're working with peripheral sub-parts, we calculate
             % position distributions as well.
             samples = (prevPosition(relevantChildren, :) - instancePositions) / halfRFSize;
             [uniqueSamples, ~, IC] = unique(samples, 'rows');
             halfPixelArr = ones(size(uniqueSamples,1),1, 'single') * halfPixel;
             
             curRange = ((instanceItr-1)*2+1):((instanceItr)*2);
             lowProbs = mvncdf(uniqueSamples + [-halfPixelArr, -halfPixelArr], mu(:,curRange), Sigma(curRange, curRange));
             leftProbs = mvncdf(uniqueSamples + [halfPixelArr, -halfPixelArr], mu(:,curRange), Sigma(curRange, curRange));
             rightProbs = mvncdf(uniqueSamples + [-halfPixelArr, halfPixelArr], mu(:,curRange), Sigma(curRange, curRange));
             highProbs = mvncdf(uniqueSamples + [halfPixelArr, halfPixelArr], mu(:,curRange), Sigma(curRange, curRange));
             pointProbs = (highProbs - (leftProbs + rightProbs)) + lowProbs;
             pointProbs(pointProbs < minProb) = minProb;
             pointProbs = pointProbs(IC);
             childrenPosActivations(:, instanceItr) = log(pointProbs);
             
             % Save samples (positions and labels).
             relevantSamples(:, ((instanceItr-1)*2+1):((instanceItr)*2)) = samples;
        end
        
%         % Finally, we choose a pdf threshold for this sub, based on
%         % the values of individual instances.
%         posActivations = pdf(obj, relevantSamples);
% 
%         % Convert to log probabilities
%         activationDenom = max(pdf(obj, obj.mu)) * 1/maxPosProb;
%         posActivations = posActivations / activationDenom;
%         posActivations = single(log(posActivations));
        
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
        newActivationArr{vocabLevelItr} = avgActivations;
        idxArr{vocabLevelItr} = idx;
        
        % Save the minimum log value for further filtering.
        vocabLevel(vocabLevelItr).minActivationLog = min(avgActivations);
    end
    
    % Finally, assign activations back.
    for vocabLevelItr = 1:numel(vocabLevel)
         newActivations(idxArr{vocabLevelItr}) = newActivationArr{vocabLevelItr};
    end
    newActivations = num2cell(newActivations);
    [graphLevel.activation] = deal(newActivations{:});
end

