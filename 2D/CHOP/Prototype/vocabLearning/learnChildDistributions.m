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
function [vocabLevel, nodeDistributionLevel] = learnChildDistributions(vocabLevel, graphLevel, prevRealLabelIds, prevPrecisePositions, prevPositions, levelItr, options)
    numberOfNodes = numel(vocabLevel);
    labelIds = [graphLevel.realLabelId];
    noiseSigma = 0.0001;
    dummySigma = 0.001;
    posDim = size(prevPrecisePositions,2);
    tempSize = getRFSize(options, levelItr);
    halfRFSize = round(tempSize(1)/2);
    maxPosProb = 0.1;
    
    % Get Receptive field size (reduced). 
    if ismembc(levelItr, options.smallRFLayers)
        rfSize = round(options.smallReceptiveFieldSize/2);
    else
        rfSize = round(options.receptiveFieldSize/2);
    end
    
    colormap('hot');
                  
    %% Distribution parameters.
    sampleCountPerCluster = 50;
    minSampleCountPerCluster = 3;
    maxClusters = 5;
    maxPointsToCluster = maxClusters * sampleCountPerCluster;
    smallSampleCount = 5;
    debugFlag = ~options.fastStatLearning && usejava('jvm');
    
    % Allocate space for distribution data.
    nodeDistributionLevel(numberOfNodes) = NodeDistribution();
    
    % Points to visualize on screen.
    pointSymbols = {'mo', 'bx', 'r*', 'kx', 'go', 'gx', 'ro', 'm*', 'mx', 'bo'};
    
    debugFolder = [options.debugFolder '/level' num2str(levelItr) '/jointStats'];
    if ~exist(debugFolder, 'dir')
         mkdir(debugFolder);
    end
    
    vocabInstances = cell(numberOfNodes,1);
    for vocabItr = 1:numberOfNodes
         vocabInstances{vocabItr} = graphLevel(labelIds == vocabItr);
    end   
    
    % Second, we go through each node, and collect statistics.
    for vocabItr = 1:numberOfNodes
        
        minX = -1; maxX = -1; minY = -1; maxY = -1; meanArr = []; scores=[]; coeff=[];
        w = warning('off', 'all');
        
        % Get data. 
        % TODO: Consider mapping here!
        instances = vocabInstances{vocabItr};
        instancePositions = cat(1, instances.precisePosition);
        instanceChildren = cat(1, instances.children);
        instanceCenterPositions = cat(1, prevPositions(instanceChildren(:,1),:));
        instanceChildrenCombinedLabels = double(prevRealLabelIds(instanceChildren)); %#ok<PFBNS>
        if size(instanceChildren,1) == 1 && size(instanceChildren,2) > 1
             instanceChildrenCombinedLabels = instanceChildrenCombinedLabels';
        end
        instanceChildrenCombinedPos = zeros(size(instanceChildren,1), size(instanceChildren,2) * posDim);
        instanceChildrenCombinedLowResPos = zeros(size(instanceChildren,1), (size(instanceChildren,2)-1) * posDim);
        
        % Find real-valued combinations. We learn a discrete distribution
        % of label combinations by simply counting the number of times each
        % combination has been seen, and then normalizing by total instance
        % count. Then, for each combination, we learn a joint distribution of
        % children positions.
        [combs, IA, IC] = unique(instanceChildrenCombinedLabels, 'rows', 'stable');
        
        % Learn mean positions of the children.
        for instanceItr = 1:size(instanceChildren,2)
             relevantChildren = instanceChildren(:, instanceItr);
             
             % If we're working with peripheral sub-parts, we calculate
             % position distributions as well.
             samples = prevPrecisePositions(relevantChildren, :) - instancePositions; %#ok<PFBNS>

             % Save samples (positions and labels).
             instanceChildrenCombinedPos(:, ((instanceItr-1)*posDim+1):((instanceItr)*posDim)) = samples;
             if instanceItr > 1
                lowResSamples = prevPositions(relevantChildren, :) - instanceCenterPositions;
                instanceChildrenCombinedLowResPos(:, ((instanceItr-2)*posDim+1):((instanceItr-1)*posDim)) = lowResSamples;
             end
        end
        
        %% Finally, we learn a joint distibution of children for every combination of real id labels. 
        % This step is crucial in having the right distributions for
        % the generative model.  For every combination, we learn a continuous label distribution.
        % Allocate space for discrete combinations, and their position distributions. 
        childrenLabelDistributions = zeros(size(combs,1), size(combs,2) + 1, 'single');
        childrenPosDistributionModes = ones(size(combs,1), 1, 'uint8');
        if numel(IA)>1
             weights = hist(IC, 1:numel(IA))';
             weights = weights / sum(weights);
        else
             weights = 1;
        end
        childrenLabelDistributions(:, 1:size(combs,2)) = combs;
        childrenLabelDistributions(:, end) = weights;

        % For every combination, learn relevant samples.
        relevantSamples = instanceChildrenCombinedPos;
        
        % Normalize relevant samples.
        relevantSamples = relevantSamples/halfRFSize;
        
        % Get low resolution samples.
        relevantLowResSamples = instanceChildrenCombinedLowResPos / rfSize;
        
        if debugFlag && ~isempty(relevantSamples)
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
             figure('Visible', 'off', 'Renderer', 'painters'), hold on;
             axis square
             subplot(1,3,1), plot(scores(:,1), scores(:,2), 'ro');
             xlim([minX, maxX]);
             ylim([minY, maxY]);
             title('Original data points')
             % Save samples.
             parsave([debugFolder '/part' num2str(vocabItr) '_data.mat'], relevantSamples, 'relevantSamples');
        end

        % If we have too many samples, it makes sense to
        % reduce number for efficiency.
        orgRelevantSamples = relevantSamples;
        if size(relevantSamples,1)>maxPointsToCluster
             idx = randperm(size(relevantSamples,1));
             idx = idx(1:maxPointsToCluster);
             orgRelevantSamples = relevantSamples;
             relevantSamples = relevantSamples(idx,:);
             if ~isempty(scores)
                scores = scores(idx,:);
             end
             relevantLowResSamples = relevantLowResSamples(idx, :);
        end
        
        minPosActivationLog = single(0);
        if ~isempty(relevantSamples)
             %% Now, we try to fit multiple gaussians to the sample points. 
             % First, we learn the number of clusters.
             noiseArr = normrnd(0, noiseSigma, size(relevantSamples));
             relevantSamples = relevantSamples + noiseArr;
             
             % Find the number of clusters, and fit a gaussian mixture
             % distribution with the given number of clusters. This is
             % a very ugly piece of code that decides numerous
             % parameters for efficiency.
%              if size(relevantSamples,1) < smallSampleCount || size(relevantSamples,1) <= size(relevantSamples,2)
%                   % We simply don't have enough data here. We switch
%                   % to simpler fitting techniques.
%                   % Learn covariance matrices / mu rows for the data.
%                   if size(relevantSamples,1) < 3
%                        covMat = eye(size(relevantSamples,2)) * dummySigma;
%                   else
%                        covMat = cov(relevantSamples);
%                        diagMatrix = eye(size(relevantSamples,2)) > 0;
%                        % Remove non-diagonal elements to simplify assumptions.
%                        covMat(~diagMatrix) = 0;
%                        covMat(diagMatrix) =  max(covMat(diagMatrix), dummySigma);
%                   end
%                   mu = mean(relevantSamples,1);
%                   
%                   % Make sure covariance matrix is valid. If not, make it
%                   % valid.
%                   try
%                        dummyProb = mvnpdf(mu, mu, covMat); %#ok<NASGU>
%                   catch %#ok<CTCH>
%                        covMat = covMat .* eye(size(covMat,1));
%                   end
%                   
%                   % If there's a very limited number of samples ( 1 or
%                   % 2), we reduce number of generated samples.
%                   % There's an intermediate number of samples,
%                   % and we can learn a proper gaussian from
%                   % these.
%                   numberOfClusters = 1;
%                   ids = ones(size(relevantSamples,1),1);
%                   obj = gmdistribution(mu, covMat);
%                   sampleCount = sampleCountPerCluster;
%              else
%                  % Here, we perform proper statistical learning.
%                  % Number of clusters are found, and clusters are
%                  % given as initializers to the gaussian fitting
%                  % process.
%                  clusteredSamples = relevantSamples;
%                  curMaxClusters = max(1, min(maxClusters, round(size(clusteredSamples,1)/minSampleCountPerCluster)));
%                  ids = mec(clusteredSamples, 'c', curMaxClusters);
%                   
%                  % Once we obtain the number of clusters, cluster points
%                  % using k-means.
%                  tempNumberOfClusters = max(ids);
%                  if tempNumberOfClusters > 1
%                       ids = kmeans(clusteredSamples, tempNumberOfClusters, 'EmptyAction', 'singleton');
%                  end
%                  numberOfClusters = max(ids);
%                   
%                  % Learn Mu+Sigma for every distribution, and initialize a
%                  % mixed distribution with the learned parameters.
%                  muArr = zeros(numberOfClusters, size(clusteredSamples,2));
%                  covArr = zeros(size(clusteredSamples,2), size(clusteredSamples,2), numberOfClusters);
%                  for clusterItr = 1:numberOfClusters
%                       muArr(clusterItr,:) = mean(clusteredSamples(ids == clusterItr,:),1);
%                       relevantClusterSamples = clusteredSamples(ids==clusterItr,:);
%                       
%                       % Obtain a covariance matrix.
%                       if size(relevantClusterSamples,1) < 3
%                             covMat = eye(size(relevantClusterSamples,2)) * dummySigma;
%                       else
%                             covMat = cov(relevantClusterSamples);
%                             diagMatrix = eye(size(relevantClusterSamples,2)) > 0 ;
%                             covMat(diagMatrix) =  max(covMat(diagMatrix), dummySigma);
%                       end
%                        
%                       % Make covMat valid.
%                       try
%                            dummyProb = mvnpdf(muArr(clusterItr,:), muArr(clusterItr,:), covMat); %#ok<NASGU>
%                       catch %#ok<CTCH>
%                            covMat = covMat .* eye(size(covMat,1));
%                       end
%                        
% %                      covMat(~diagMatrix) = 0;
%                       
%                       % Finally, save covariance matrix.
%                       covArr(:,:,clusterItr) = covMat;
%                  end
%                   sampleCount = sampleCountPerCluster * numberOfClusters;
%                  
%                   % Calculate weights.
%                   p = hist(ids, unique(ids));
%                   p = p/sum(p);
%                   
%                   % Obtain the final gaussian mixture.
%  %                 obj = gmdistribution(muArr, covArr, p);
%                   
%                   % Finally, fit a multi=dimensional gaussian given the
%                   % parameters.
%                   options = statset('MaxIter', 100);
%                   obj = gmdistribution.fit(clusteredSamples, numberOfClusters, 'CovType', 'diagonal', 'Regularize', dummySigma, 'Start', ids, 'Options', options);
%              end

             % Simplified distribution learning.
             if size(relevantSamples,1) > size(relevantSamples,2)
                 try
                    obj = fitgmdist(relevantSamples,1, 'RegularizationValue', dummySigma, 'Replicates', 3);
                    dummyActivation = pdf(obj, relevantSamples(1,:)); %#ok<NASGU>
                 catch
                      try 
                          obj = gmdistribution.fit(relevantSamples,1, 'Regularize', dummySigma, 'Replicates', 3);
                          dummyActivation = pdf(obj, relevantSamples(1,:)); %#ok<NASGU>
                      catch
                          mu = mean(relevantSamples,1);
                          covMat = cov(relevantSamples);
                          diagMatrix = eye(size(relevantSamples,2)) > 0 ;
                          covMat(diagMatrix) =  max(covMat(diagMatrix), dummySigma);
                          covMat = eye(size(covMat));
                          obj = gmdistribution(mu, covMat);
                      end
                 end
             else
                 if size(relevantSamples,1) == 1
                     covMat = eye(size(relevantSamples,2)) * dummySigma;
                     mu = relevantSamples;
                 else
                     covMat = cov(relevantSamples);
                     diagMatrix = eye(size(relevantSamples,2)) > 0;
                     % Remove non-diagonal elements to simplify assumptions.
                     covMat(~diagMatrix) = 0;
                     covMat(diagMatrix) =  max(covMat(diagMatrix), dummySigma);
                     mu = mean(relevantSamples,1);
                     [~, muIdx] = min(pdist2(mu, relevantSamples));
                     mu = relevantSamples(muIdx,:);
                 end
                 obj = gmdistribution(mu, covMat);
             end
             sampleCount = sampleCountPerCluster;

             % Calculate maximum possible position activation.
             activationDenom = max(pdf(obj, obj.mu)) * 1/maxPosProb;
             
             % Generate random samples for density estimation.
             y = random(obj, sampleCount);
             yActivations = pdf(obj, y);
             yActivations = yActivations / activationDenom;
             yActivations = yActivations/max(yActivations);
             
%              % Finally, we choose a pdf threshold for this sub, based on
%              % the values of individual instances.
%              posActivations = pdf(obj, orgRelevantSamples);
%              % Convert activations into pseudo-probabilities.
%              posActivations = posActivations / activationDenom;
%              
%              % Finally, take log of the minimum value and save it.
%              minPosActivationLog = single(log(min(posActivations)));
             
             % Recalculate activations for the small set of samples.
             posActivations = pdf(obj, relevantSamples);
             % Convert activations into pseudo-probabilities.
             posActivations = posActivations / activationDenom;
             
             % Normalize activations for visualization.
             normPosActivations = posActivations/max(posActivations);
             
             if debugFlag && ~isempty(relevantSamples)
                  % Visualize clusters.
%                   for clusterItr = 1:numberOfClusters
%                         plot(scores(ids == clusterItr,1), scores(ids == clusterItr,2), pointSymbols{clusterItr}); %#ok<PFBNS>
%                   end
                  subplot(1,3,2), scatter(scores(:,1), scores(:,2), [], normPosActivations, 'filled');
                  
                  xlim([minX, maxX]);
                  ylim([minY, maxY]);
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
                  title('Uni-modal distribution') 

                  axesHandles = findobj(get(gcf,'Children'), 'flat','Type','axes');
                  axis(axesHandles,'square');
                  
                  % Stop drawing, move on.
                  hold off;
                  hardcopy(gcf, [debugFolder '/part' num2str(vocabItr) '.eps'], '-depsc2');
        %          saveas(gcf, [debugFolder '/part' num2str(vocabItr) '.png']);
                  close(gcf);
                  
                  %% Now, we will try to visualize each sub-part's distribution separately.
                  figure('Visible', 'off', 'Renderer', 'painters'), hold on;
                  axis square
                  numberOfSubParts = (size(relevantSamples,2)/2);
                  for subPartItr = 1:numberOfSubParts
                      tempSamples = relevantSamples(:, ((subPartItr-1)*2+1):subPartItr*2);
                      subplot(2, numberOfSubParts, subPartItr), scatter(tempSamples(:,2), tempSamples(:,1), [], normPosActivations, 'filled');
                      xlim([-1, 1]);
                      ylim([-1, 1]);
                      title(['Child ' num2str(subPartItr) ' (Data)']); 
                  end
                  
                  %% Visualize each sub-part of generated data this time.
                  for subPartItr = 1:numberOfSubParts
                      tempSamples = y(:, ((subPartItr-1)*2+1):subPartItr*2);
                      subplot(2, numberOfSubParts, numberOfSubParts + subPartItr), scatter(tempSamples(:,2), tempSamples(:,1), [], yActivations, 'filled');
                      xlim([-1, 1]);
                      ylim([-1, 1]);
                      title(['Child ' num2str(subPartItr) ' (Generated)']); 
                  end
                  
                  % Make the subplots square.
                  axesHandles = findobj(get(gcf,'Children'), 'flat','Type','axes');
                  axis(axesHandles,'square');
                  hold off;
                  
                  % Save figure.
                  hardcopy(gcf,[debugFolder '/part' num2str(vocabItr) '_subparts.eps'],'-depsc2');
     %             saveas(gcf, [debugFolder '/part' num2str(vocabItr) '_subparts.png']);
                  close(gcf);
             end

             % Here, we assign each label combination to the closest gaussian.
             for combItr = 1:size(combs,1)
                  combSamples = orgRelevantSamples(IC==combItr, :);
                  meanCombSamples = mean(combSamples,1);

                  % Choose the closest mode.
                  mixtureMus = obj.mu;
                  mixtureSigmas = obj.Sigma;
                  weights = obj.PComponents;
                  pdfVals = weights;
                  for mixItr = 1:size(mixtureMus,1)
                       try
                            sigmaMatrix = squeeze(mixtureSigmas(:,:,mixItr));
                       catch %#ok<CTCH>
                            sigmaMatrix = eye(numel(mixtureSigmas)) * mixtureSigmas(1);
                       end
%                      try
                       pdfVals(mixItr) = mvnpdf(meanCombSamples, mixtureMus(mixItr,:), sigmaMatrix) * pdfVals(mixItr);
%                        catch %#ok<CTCH>
%                             sigmaMatrix = sigmaMatrix .* eye(size(sigmaMatrix,1)) * noiseSigma;
%                             pdfVals(mixItr) = mvnpdf(meanCombSamples, mixtureMus(mixItr,:), sigmaMatrix) * pdfVals(mixItr);
%                        end
                  end

                  % Select most likely mode. If probability is zero, select the
                  % the mode with the most weight.
                  [val, idx] = max(pdfVals);
                  if val == 0
                       [~, idx] = max(weights);
                  end

                  % Assign mode.
                  childrenPosDistributionModes(combItr) = idx;
             end
        else
             obj = gmdistribution([], []);
        end
        
        % Assign distributions.
        nodeDistributionLevel(vocabItr).childrenLabelDistributions = childrenLabelDistributions;
        nodeDistributionLevel(vocabItr).childrenPosDistributions = obj;
 %       nodeDistributionLevel(vocabItr).minPosActivationLog = minPosActivationLog;
        nodeDistributionLevel(vocabItr).childrenPosDistributionModes = childrenPosDistributionModes;
        
        % Re-open warnings.
        warning(w);
    end
    
    colormap default;
    clearvars -except vocabLevel nodeDistributionLevel
end

function parsave(fname, data, dataName) %#ok<INUSL>
     eval([dataName '=data;']);
     save(fname, 'dataName')
end