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
function [vocabLevel] = learnChildDistributions(vocabLevel, graphLevel, previousLevel, levelItr, options)
    numberOfNodes = numel(vocabLevel);
    labelIds = [graphLevel.labelId];
    prevRealLabelIds = [previousLevel.realLabelId]';
    prevPosition = double(cat(1, previousLevel.precisePosition));
    noiseSigma = 0.0001;
    dummySigma = 0.001;
    posDim = size(prevPosition,2);
    sampleCountPerCluster = 100;
    maxPointsToCluster = 500;
    smallSampleCount = 10;
    maxClusters = 10;
    minPoints = maxClusters + 5;
    debugFlag = ~options.fastStatLearning;
    
    % Close warnings.
    
    pointSymbols = {'mo', 'bx', 'r*', 'kx', 'go', 'b*', 'ro', 'gx', 'm*', 'k*'};
    
    debugFolder = [options.debugFolder '/level' num2str(levelItr) '/jointStats'];
    if ~exist(debugFolder, 'dir')
         mkdir(debugFolder);
    end
    
    vocabInstances = cell(numberOfNodes,1);
    for vocabItr = 1:numberOfNodes
         vocabInstances{vocabItr} = graphLevel(labelIds == vocabItr);
    end                       
    
    % Second, we go through each node, and collect statistics.
    parfor vocabItr = 1:numberOfNodes
        minX = -1; maxX = -1; minY = -1; maxY = -1; meanArr = []; scores=[]; coeff=[];
        w = warning('off', 'all');
        children = vocabLevel(vocabItr).children;
        
        % Get data. 
        % TODO: Consider mapping here!
        instances = vocabInstances{vocabItr};
        instanceChildren = cat(1, instances.children);
        instanceChildrenCombinedLabels = double(prevRealLabelIds(instanceChildren)); %#ok<PFBNS>
        if size(instanceChildren,1) == 1 && size(instanceChildren,2) > 1
             instanceChildrenCombinedLabels = instanceChildrenCombinedLabels';
        end
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
                  samples = prevPosition(relevantChildren, :) - prevPosition(instanceChildren(:, 1), :); %#ok<PFBNS>
                  
                  % Save samples (positions and labels).
                  instanceChildrenCombinedPos(:, ((instanceItr-2)*posDim+1):((instanceItr-1)*posDim)) = samples;
             end
        end
        
       
        %% Finally, we learn a joint distibution of children for every combination of real id labels. 
        % This step is crucial in having the right distributions for
        % the generative model.  For every combination, we learn a continuous label distribution.
        % Allocate space for discrete combinations, and their position distributions. 
        childrenLabelDistributions = zeros(size(combs,1), size(combs,2) + 1, 'single');
        childrenPosDistributions = cell(1);
        if numel(IA)>1
             weights = hist(IC, 1:numel(IA))';
             weights = weights / sum(weights);
        else
             weights = 1;
        end
        childrenLabelDistributions(:, 1:size(combs,2)) = combs;
        childrenLabelDistributions(:, end) = weights;

        % For every combination, learn relevant samples.
        if ~isempty(instanceChildrenCombinedPos)
             relevantSamples = instanceChildrenCombinedPos;
             if debugFlag
                  if size(relevantSamples,2) == 2
                       coeff = eye(2);
                       scores = relevantSamples;
                  else
                       [coeff, scores] = princomp(relevantSamples);
                  end
                  meanArr = mean(relevantSamples,1);

                  % Find min/max values for showing sample points on screen.
                  minX = min(scores(:,1)) - 1;
                  maxX = max(scores(:,1)) + 1;
                  minY = min(scores(:,2)) - 1;
                  maxY = max(scores(:,2)) + 1;

                  %% Get principle components and visualize samples in 2d space.
                  figure('Visible', 'off'), hold on;
                  axis square
                  subplot(1,3,1), plot(scores(:,1), scores(:,2), 'ro');
                  xlim([minX, maxX]);
                  ylim([minY, maxY]);
                  title('Original data points')
                  % Save samples.
                  parsave([debugFolder '/part' num2str(vocabItr) '.mat'], relevantSamples, 'relevantSamples');
             end

             %% Now, we try to fit multiple gaussians to the sample points. 
             % First, we learn the number of clusters.
             noiseArr = normrnd(0, noiseSigma, size(relevantSamples));
             relevantSamples = relevantSamples + noiseArr;

             % Find the number of clusters, and fit a gaussian mixture
             % distribution with the given number of clusters. This is
             % a very ugly piece of code that decides numerous
             % parameters for efficiency.
             regularizeTerm = 1e-5;
             bigRegularizeTerm = 0.1;
             if size(relevantSamples,1) < minPoints || size(relevantSamples,1) <= size(relevantSamples,2)
                  % We simply don't have enough data here. We switch
                  % to simpler fitting techniques.
                  % Learn covariance matrices / mu rows for the data.
                  if size(relevantSamples,1) == 1
                       covMat = dummySigma;
                       mu = relevantSamples;
                  else
                       if size(relevantSamples,1) <= size(relevantSamples,2)
                             covMat = cov(relevantSamples) + bigRegularizeTerm;
                       else
                             covMat = cov(relevantSamples) + regularizeTerm;
                       end
                       mu = mean(relevantSamples,1);
                  end

                  % If there's a very limited number of samples ( 1 or
                  % 2), we reduce number of generated samples.
                  if size(relevantSamples,1) < 3
                       numberOfClusters = size(relevantSamples,1);
                       ids = 1:numberOfClusters;
                       obj = gmdistribution(relevantSamples, repmat(dummySigma, 1, size(relevantSamples,2)));
                       if size(relevantSamples,1) == 1
                            sampleCount = 1;
                       else
                            sampleCount = smallSampleCount;
                       end
                  else
                       % There's an intermediate number of samples,
                       % and we can learn a proper gaussian from
                       % these.
                       numberOfClusters = 1;
                       ids = ones(size(relevantSamples,1),1);
                       obj = gmdistribution(mu, covMat);
                       sampleCount = sampleCountPerCluster;
                  end
             else

                  % If we have too many samples, it makes sense to
                  % reduce number for efficiency.
                  if size(relevantSamples,1)>maxPointsToCluster
                       idx = randperm(size(relevantSamples,1));
                       idx = idx(1:maxPointsToCluster);
                       relevantSamples = relevantSamples(idx,:);
                  end

                  % Here, we perform proper statistical learning.
                  % Number of clusters are found, and clusters are
                  % given as initializers to the gaussian fitting
                  % process.
                  ids = mec(relevantSamples, 'c', maxClusters);
                  numberOfClusters = max(ids);
                  sampleCount = sampleCountPerCluster * numberOfClusters;
                  options = statset('MaxIter', 5);
                  obj = gmdistribution.fit(relevantSamples, numberOfClusters, 'Regularize', regularizeTerm, 'Start', ids, 'Options', options);
             end


             % Save distributions.
             childrenPosDistributions = {obj};

             if debugFlag
                  % Generate random samples for density estimation.
                  y = random(obj, sampleCount);
                  
                  % Visualize clusters.
                  subplot(1,3,2), hold on 
                  for clusterItr = 1:numberOfClusters
                        plot(scores(ids == clusterItr,1), scores(ids == clusterItr,2), pointSymbols{clusterItr}); %#ok<PFBNS>
                  end
                  xlim([minX, maxX]);
                  ylim([minY, maxY]);
                  hold off
                  title('Seed clusters for Gaussian Mixtures');

                  % Calculate new scores for visualization.
                  if size(relevantSamples,2) == 2
                       newScores = y;
                  else
                       newScores = (y - repmat(meanArr, size(y,1), 1)) * coeff;
                       newScores = newScores(:,1:2);
                  end
                  subplot(1,3,3), plot(newScores(:,1), newScores(:,2), 'mo');
                  xlim([minX, maxX]);
                  ylim([minY, maxY]);
                  title('Multi-modal (max 5) distribution') 

                  % Stop drawing, move on.
                  hold off
                  saveas(gcf, [debugFolder '/part' num2str(vocabItr) '.png']);
                  close(gcf);
             end
        end
        
        vocabLevel(vocabItr).childrenLabelDistributions = childrenLabelDistributions;
        vocabLevel(vocabItr).childrenPosDistributions{1} = childrenPosDistributions;
        vocabLevel(vocabItr).realChildren = children;
        
        % Re-open warnings.
        warning(w);
    end
    clearvars -except vocabLevel
end

function parsave(fname, data, dataName) %#ok<INUSL>
     eval([dataName '=data;']);
     save(fname, 'dataName')
end