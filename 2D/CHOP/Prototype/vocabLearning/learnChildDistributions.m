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
function [vocabLevel] = learnChildDistributions(vocabLevel, graphLevel, previousLevel, previousModes, modeProbs, levelItr, options)
    numberOfNodes = numel(vocabLevel);
    labelIds = [graphLevel.labelId];
    prevRealLabelIds = [previousLevel.realLabelId];
    prevPosition = double(cat(1, previousLevel.position));
    dummySigma = 0.001;
    posDim = size(prevPosition,2);
    generatedSampleCount = 200;
    halfSize = (options.receptiveFieldSize/2) + 3;
    midPoint = floor(options.receptiveFieldSize/2) + 1;
    midPoint = [midPoint, midPoint];
    maxFitTrials = 10;
    
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
             
             % If we're working with peripheral sub-parts, we calculate
             % position distributions as well.
             if instanceItr > 1
                  samples = prevPosition(relevantChildren, :) - prevPosition(instanceChildren(:, 1), :);
                  
                  % Save samples (positions and labels).
                  instanceChildrenCombinedPos(:, ((instanceItr-2)*posDim+1):((instanceItr-1)*posDim)) = samples;

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
        if size(instanceChildrenCombinedPos,1) > 5 && numel(children)>1
             for combItr = 1:numel(IA)
                  relevantSamples = instanceChildrenCombinedPos(IC==combItr,:);
                  [coeff, score] = princomp(relevantSamples);
                  meanArr = mean(relevantSamples,1);
                  
                  % Save samples.
                  save([debugFolder '/part' num2str(vocabItr) '_' mat2str(combs(combItr,:)) '.mat'], 'relevantSamples');
                  
                  %% First,, draw borders of the area that can be represented.
                  % First, we obtain all possible positions
                  % for every child. 
                  positionArr = cell(size(instanceChildren,2)-1, 1);
                  for itr = 1:(size(instanceChildren,2) - 1)
                       curEdge = vocabLevel(vocabItr).adjInfo(itr,:);
                       curEdge(1:2) = vocabLevel(vocabItr).children(curEdge(1:2));
                       relevantMode = ismember(previousModes(:,1:3), curEdge(1:3), 'rows');
                       positions = find(squeeze(modeProbs(relevantMode, :, :)));
                       positionArr{itr} = positions;
                  end
                  allCombs = allcomb(positionArr{:});
                  
                  % Create samples corresponding to all combinations.
                  allCombSamples = zeros(size(allCombs,1), size(allCombs,2) * posDim);
                  for itr = 1:size(allCombs,2)
                       [allCombSamplesX, allCombSamplesY] = ind2sub([options.receptiveFieldSize, options.receptiveFieldSize], allCombs(:,itr));
                       allCombSamples(:, ((itr-1) * posDim + 1):(itr * posDim)) = [allCombSamplesX, allCombSamplesY] - repmat(midPoint, size(allCombs,1), 1);
                  end
                  
                  % Calculate scores for the combinations.
                  newScores = (allCombSamples - repmat(meanArr, size(allCombSamples,1), 1)) * coeff;
                  newScores = newScores(:,1:2);
                  try
                       K = convhull(newScores(:,1), newScores(:,2));
                  catch %#ok<CTCH>
                       % Colinear points.
                       K = ones(size(newScores,1),1)>0;
                  end
                  convHullPoints = newScores(K,:);
                  figure, hold on
                  subplot(2,3,1), plot(convHullPoints(:,1), convHullPoints(:,2), 'r-', newScores(:,1), newScores(:,2), 'b*');
                  minX = min(convHullPoints(:,1)) - 1;
                  maxX = max(convHullPoints(:,1)) + 1;
                  minY = min(convHullPoints(:,2)) - 1;
                  maxY = max(convHullPoints(:,2)) + 1;
                  xlim([minX, maxX]);
                  ylim([minY, maxY]);
                  title('Area covered by all possible examples (Expressivity)') 

                  %% Get principle components and visualize samples in 2d space.
                  axis square
                  subplot(2,3,2), plot(convHullPoints(:,1), convHullPoints(:,2), 'r-', score(:,1), score(:,2), 'ro');
                  xlim([minX, maxX]);
                  ylim([minY, maxY]);
                  title('Original data points')  
                  
                  %% As a second option, we sample from individual distributions.
                  newSamples = zeros(generatedSampleCount, size(relevantSamples, 2));
                  for instanceItr = 2:size(instanceChildren,2)
                        newSamples(:,((instanceItr-2)*posDim+1):((instanceItr-1)*posDim)) = ...
                             mvnrnd(meanChildrenPos(instanceItr,:), ...
                             [meanChildrenCov(instanceItr,1:2); meanChildrenCov(instanceItr,3:4)], generatedSampleCount);
                  end
                  newScores = (newSamples - repmat(meanArr, size(newSamples,1), 1)) * coeff;
                  newScores = newScores(:,1:2);
                  subplot(2,3,4), plot(convHullPoints(:,1), convHullPoints(:,2), 'r-', newScores(:,1), newScores(:,2), 'bo');
                  xlim([minX, maxX]);
                  ylim([minY, maxY]);
                  title('Samples (Individual distributions)')  
                  
                  %% The third option is to fit a gaussian to the joint distribution and sample from that.
                  % First, we calculate a standard deviation for every
                  % column. If in any column it's zero, we add gaussian
                  % noise to the values.
                  devs = std(relevantSamples, 0, 1);
                  if ismember(0, devs)
                       noiseArr = normrnd(0, dummySigma, size(relevantSamples));
                       relevantSamples = relevantSamples + noiseArr;
                  end
                  
                  % Fit distribution and generate samples.
                  try
                       obj = gmdistribution.fit(relevantSamples,1);
                       y = random(obj, generatedSampleCount);
                  catch %#ok<CTCH>
                       mu = mean(relevantSamples,1);
                       covMat = cov(relevantSamples);
                       y = mvnrnd(mu, covMat, generatedSampleCount);
                  end
                  newScores = (y - repmat(meanArr, size(y,1), 1)) * coeff;
                  newScores = newScores(:,1:2);
                  subplot(2,3,5), plot(convHullPoints(:,1), convHullPoints(:,2), 'r-', newScores(:,1), newScores(:,2), 'bo');
                  xlim([minX, maxX]);
                  ylim([minY, maxY]);
                  title('Samples (Joint distribution)') 
                  
                  %% Now, we try to fit multiple gaussians to the sample points. 
                  % First, we learn the number of clusters.
                  noiseArr = normrnd(0, dummySigma, size(relevantSamples));
                  relevantSamples = relevantSamples + noiseArr;
                  
                  ids = mec(relevantSamples, 'c', 3);
                  numberOfClusters = max(ids);
                  
                  obj = gmdistribution.fit(relevantSamples, numberOfClusters, 'Regularize', 1e-10);
                  y = random(obj, generatedSampleCount);
                  
%                  if numberOfClusters > 1                    
%                        fitFlag = false;
%                        trialCounter = 0;
%                        while ~fitFlag && trialCounter < maxFitTrials
%                             fitFlag = true;
%                             try
%                                  obj = gmdistribution.fit(relevantSamples,numberOfClusters);
%                                  y = random(obj, generatedSampleCount);
%                             catch 
%                                  fitFlag = false;
%                             end
%                        end
                        
%                        if ~fitFlag
%                             fitFlag = 0;
%                             sampleCounts = round((hist(ids, unique(ids)) / numel(ids)) * generatedSampleCount);
%                             muArr = cell(numberOfClusters,1);
%                             covArr = cell(numberOfClusters,1);
%                             for clusterItr = 1:numberOfClusters
%                                  idx = ids == clusterItr;
%                                  muArr{clusterItr} = mean(relevantSamples(ids == clusterItr,:),1);
%                                  if nnz(idx) > 1
%                                       covArr{clusterItr} = cov(relevantSamples(ids == clusterItr,:));
%                                  else
%                                       dummyCov = zeros(size(relevantSamples,2));
%                                       dummyCov(1:(size(relevantSamples,2)+1):(size(relevantSamples,2)^2)) = dummySigma;
%                                       covArr{clusterItr} = dummyCov;
%                                  end
%                             end
%                             
%                             % Generate samples from individual
%                             % distributions.
%                             y = zeros(sum(sampleCounts), size(relevantSamples,2));
%                             offset = 1;
%                             for clusterItr = 1:numberOfClusters
%                                  y(offset:(offset+sampleCounts(clusterItr)-1), :) = mvnrnd(muArr{clusterItr}, covArr{clusterItr}, sampleCounts(clusterItr));
%                                  offset = offset + sampleCounts(clusterItr);
%                             end
%                        end
%                  end
                  newScores = (y - repmat(meanArr, size(y,1), 1)) * coeff;
                  newScores = newScores(:,1:2);
                  subplot(2,3,6), plot(convHullPoints(:,1), convHullPoints(:,2), 'r-', newScores(:,1), newScores(:,2), 'bo');
                  xlim([minX, maxX]);
                  ylim([minY, maxY]);
                  title('Multi-modal (max 3) distribution') 
                  
                  
                  %% T-SNE results!
 %                 mappedX = tsne(relevantSamples, [], 2, size(relevantSamples,2), 5);
 %                 subplot(2,3,3), plot(convHullPoints(:,1), convHullPoints(:,2), 'r-', mappedX(:,1), mappedX(:,2), 'ko');
  %                title('t-SNE neighbor embedding') 
                  % Stop drawing, move on.
                  hold off
                  
                  % Find principal components and visualize them.
                  saveas(gcf, [debugFolder '/part' num2str(vocabItr) '_' mat2str(combs(combItr,:)) '.png']);
                  close(gcf);
             end
        end
        
        vocabLevel(vocabItr).childrenPosMean = meanChildrenPos;
        vocabLevel(vocabItr).childrenPosCov = meanChildrenCov;
        vocabLevel(vocabItr).realChildren = children;
    end
end