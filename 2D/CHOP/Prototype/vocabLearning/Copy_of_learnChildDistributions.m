%> Name: learnChildDistributions
%>
%> Description: This function learns label/position distributions of the sub-parts
%> of every part in vocabLevel, using the data in the object graphs.
%> 
%> @param vocabLevel Current vocabulary level.
%> @param graphLevel Current graph level.
%> @param previousLevel Previous graph level.
%> 
%> @retval vocabLevel The vocabulary level updated with sub-part
%> distributions.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 07.07.2015
function [vocabLevel] = learnChildDistributions(vocabLevel, graphLevel, previousLevel, previousModes, levelItr, options)
    numberOfNodes = numel(vocabLevel);
    labelIds = [graphLevel.labelId];
    curPrecisePos = cat(1, graphLevel.precisePosition);
    prevRealLabelIds = [previousLevel.realLabelId];
    prevPrecisePos = cat(1, previousLevel.precisePosition);
    dummySigma = 0.1;
    posDim = size(prevPrecisePos,2);
    generatedSampleCount = 200;
    halfSize = (options.receptiveFieldSize/2) + 1;
    
    debugFolder = [options.debugFolder '/level' num2str(levelItr) '/jointStats'];
    if ~exist(debugFolder, 'dir')
         mkdir(debugFolder);
    end
    
    % Second, we go through each node, and collect statistics.
    for vocabItr = 1:numberOfNodes
        children = vocabLevel(vocabItr).children;
        curPart = vocabLevel(vocabItr);
        meanChildrenPos = zeros(numel(children),2, 'single');
        meanChildrenCov = zeros(numel(children),4,'single');
        
        % Get data. 
        % TODO: Consider mapping here!
        relevantIdx = labelIds == vocabItr;
        instances = graphLevel(relevantIdx);
        instanceChildren = cat(1, instances.children);
        instanceChildrenCombinedLabels = double(prevRealLabelIds(instanceChildren));
        instanceChildrenCombinedPos = zeros(size(instanceChildren,1), (size(instanceChildren,2)-1) * posDim);
        realIdInstanceChildren  = prevRealLabelIds(instanceChildren);
        
        % Find real-valued combinations. We learn a discrete distribution
        % of label combinations by simply counting the number of times each
        % combination has been seen, and then normalizing by total instance
        % count. Then, for each combination, we learn a joint distribution of
        % children positions.
        [combs, IA, IC] = unique(instanceChildrenCombinedLabels, 'rows', 'stable');
       
        % Learn mean positions of the children.
        for instanceItr = 1:size(instanceChildren,2)
             relevantChildren = instanceChildren(:, instanceItr);
             childLabel = mode(double(prevRealLabelIds(relevantChildren)));
             children(instanceItr) = childLabel;
             samples = prevPrecisePos(relevantChildren, :) - curPrecisePos(relevantIdx, :);
             
             % Save samples (positions and labels).
             instanceChildrenCombinedPos(:, ((instanceItr-1)*posDim+1):(instanceItr*posDim)) = samples;
             
             % Calculate mean of the distribution (Gaussian assumption).
             meanChildrenPos(instanceItr,:) = mean(samples,1);
             
             % Calculate the covariance of the samples. 
             if size(samples,1) > 1
                 covMat = cov(samples);
                 meanChildrenCov(instanceItr,1:2) = covMat(1,:);
                 meanChildrenCov(instanceItr,3:4) = covMat(2,:);
                 if nnz(meanChildrenCov(instanceItr,1:4)) < 2
                    meanChildrenCov(instanceItr,[1, 4]) = max(meanChildrenCov(instanceItr,[1, 4]), [dummySigma, dummySigma]);
                 end
             else
                 meanChildrenCov(instanceItr,[1, 4]) = dummySigma;
             end
        end
        
        %% Finally, we learn a joint distibution of children for every combination of real id labels. 
        % This step is crucial in having the right distributions for
        % the generative model. The following actions are essentially
        % for debugging. First, we generate samples for the most
        % generic case, which is what the vocabulary parts are capable
        % of representing. We would like to show that the grammar is
        % overly-expressive. 

        % First, we start by the basic model: What's capable of
        % activating this part? We generate random samples from the
        % part description. 
        
        % For every combination, we learn a continuous label distribution.
        if size(instanceChildrenCombinedPos,1) > 5
             for combItr = 1:numel(IA)
                  relevantSamples = instanceChildrenCombinedPos(IC==combItr,:);

                  %% Get principle components and visualize samples in 2d space.
                  [coeff, score] = princomp(relevantSamples);
                  figure,
                  subplot(1,2,1), plot(score(:,1), score(:,2), 'ro');
                  xlim([-halfSize, halfSize]);
                  ylim([-halfSize, halfSize]);
                  
                  %% As a second option, we sample from individual distributions.
                  newSamples = zeros(size(relevantSamples));
                  for instanceItr = 1:size(instanceChildren,2)
                        newSamples(:,((instanceItr-1)*posDim+1):(instanceItr*posDim)) = ...
                             mvnrnd(meanChildrenPos(instanceItr,:), ...
                             [meanChildrenCov(instanceItr,1:2); meanChildrenCov(instanceItr,3:4)], size(relevantSamples,1));
                  end
                  
                  % Find principal components and visualize them.
                  newScores = (newSamples - repmat(mean(newSamples,1), size(newSamples,1), 1)) * coeff;
                  newScores = newScores(:,1:2);
                  subplot(1,2,2), plot(newScores(:,1), newScores(:,2), 'bo');
                  xlim([-halfSize, halfSize]);
                  ylim([-halfSize, halfSize]);
                  saveas(gcf, [debugFolder '/part' num2str(vocabItr) '_' mat2str(combs(combItr,:)) '.png']);
                  close(gcf);
             end
        end
        
        vocabLevel(vocabItr).childrenPosMean = meanChildrenPos;
        vocabLevel(vocabItr).childrenPosCov = meanChildrenCov;
        vocabLevel(vocabItr).realChildren = children;
    end
end