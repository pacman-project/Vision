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
function [vocabLevel, nodeDistributionLevel] = learnChildDistributions(vocabLevel, nodeDistributionLevel, graphLevel, prevRealLabelIds, prevPrecisePositions, prevPositions, levelItr, options)
    numberOfNodes = numel(vocabLevel);
    labelIds = [graphLevel.realLabelId];
    posDim = size(prevPrecisePositions,2);
    tempSize = getRFSize(options, levelItr);
    realRFSize = tempSize(1);
    centerRF = floor(tempSize(1)/2)+1;
    lowRFSize = options.receptiveFieldSizes(levelItr-1);
    minProb = exp(-20);
    distributions = {nodeDistributionLevel.childrenPosDistributions};
    
    colormap('hot');

    %% Distribution parameters.
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
        w = warning('off', 'all');
        
        % Get data. 
        % TODO: Consider mapping here!
        instances = vocabInstances{vocabItr};
        instancePositions = cat(1, instances.precisePosition);
        instanceChildren = cat(1, instances.children);
        instanceChildrenCombinedLabels = double(prevRealLabelIds(instanceChildren));
        if size(instanceChildren,1) == 1 && size(instanceChildren,2) > 1
             instanceChildrenCombinedLabels = instanceChildrenCombinedLabels';
        end
        instanceChildrenCombinedPos = zeros(size(instanceChildren,1), size(instanceChildren,2) * posDim);
        
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
             samples = prevPrecisePositions(relevantChildren, :) - instancePositions;

             % Save samples (positions and labels).
             instanceChildrenCombinedPos(:, ((instanceItr-1)*posDim+1):((instanceItr)*posDim)) = samples;
        end
        
        %% Finally, we learn a joint distibution of children for every combination of real id labels. 
        % This step is crucial in having the right distributions for
        % the generative model.  For every combination, we learn a continuous label distribution.
        % Allocate space for discrete combinations, and their position distributions. 
        childrenLabelDistributions = zeros(size(combs,1), size(combs,2) + 1, 'single');
        if numel(IA)>1
             weights = hist(IC, 1:numel(IA))';
             weights = weights / sum(weights);
        else
             weights = 1;
        end
        childrenLabelDistributions(:, 1:size(combs,2)) = combs;
        childrenLabelDistributions(:, end) = weights;

        %% After we have learned the distributions, now let us calculate precise sub-children 
        % probabilities. This will be used when we are calculating
        % probabilities for each sub-location.
        % Start by creating real-sized receptive fields.
        childrenProbs = cell(size(instanceChildren,2),1);
        childImgs = zeros(realRFSize, size(instanceChildren,2) * realRFSize, 3, 'uint8');
        
        % Next, we enumerate all possible (Really!) points.
        oneSideVals = (-(centerRF-1)):1:(realRFSize - centerRF);
        minRF = (-(centerRF-1));
        maxRF = (realRFSize - centerRF);
        combs = allcomb(oneSideVals, oneSideVals);
        
        % We calculate pdf values. For locations with significiant pdf
        % values, we will evaluate cdf (since it is too costly!).
        obj = distributions{vocabItr};
        mu = obj.mu;
        Sigma = obj.Sigma;
        for childItr = 1:size(instanceChildren,2)
             % Get relevant distribution, and calculate pdf values for
             % early filtering.
             curMu = mu(:,((childItr-1)*2+1):(childItr*2));
             curSigma = Sigma(((childItr-1)*2+1):(childItr*2),((childItr-1)*2+1):(childItr*2));
             try
                  combPDFVals = mvnpdf(combs, curMu, curSigma);
             catch
                  combPDFVals = mvnpdf(combs, curMu, nearestSPD(curSigma));
             end
             pdfMat = zeros(realRFSize, realRFSize, 'double');
             pdfMat(combs(:,1)+centerRF + (combs(:,2)+centerRF-1)*realRFSize) = combPDFVals;
             
             % Determine the value of minimum pdf to check.
             relevantPos = unique(instanceChildrenCombinedPos(:,((childItr-1)*2 + 1):(childItr * 2)), 'rows');
             relevantPos = relevantPos(relevantPos(:,1) >= minRF & relevantPos(:,1) <= maxRF & ...
                  relevantPos(:,2) >= minRF & relevantPos(:,2) <= maxRF, :);
             relevantInd = sub2ind([realRFSize, realRFSize], relevantPos(:,1) + centerRF, relevantPos(:,2) + centerRF);
             sampleImg = zeros(realRFSize) > 0;
             sampleImg(relevantInd) = 1;
             dummyArr = sparse(realRFSize,realRFSize);
             
             % Based on the number of instances, we add a buffer zone in
             % probability distribution to add flexibility.
             if levelItr > 2
                  bufferWidth = max(1, round(realRFSize /( 2*(lowRFSize * sqrt(size(relevantPos,1))))));
                  sampleImg = imdilate(sampleImg, strel('disk', bufferWidth));
             end
             
             % Finally, we calculate a convex hull and fill it in.
             sampleImg = bwconvhull(sampleImg);
             
             % Find points and save them.
             if isempty(relevantPos)
                  childrenProbs{childItr} = dummyArr;
                  continue;
             end
             
             pdfMat(~sampleImg) = 0;
             % We calculate pseudo-probabilities here. Simply sum all
             % valid pdf values, and then calculate probability by dividing them
             % by their sum.
             dummyArr = pdfMat / sum(sum(pdfMat));
             dummyArr(dummyArr < minProb & dummyArr > 0) = minProb;
             
             % Finally, assign the calculated probabilities to the correct
             % positions.
             childrenProbs{childItr} = sparse(dummyArr);
             
             % Create a visualization.
             maxVal = max(max(dummyArr));
             dummyArr = dummyArr / maxVal;
             dummyArr([1, end], :) = 1;
             dummyArr(:, [1, end]) = 1;
 %            relevantPos = relevantPos + centerRF;
 %            ind = sub2ind([realRFSize, realRFSize], relevantPos(:,1), relevantPos(:,2));
 %            dummyArr(ind) = 1;
             dummyArr(dummyArr>0) = max(1/255, dummyArr(dummyArr>0));
             img = label2rgb(uint8(255 * dummyArr), 'jet', 'k');
             childImgs(:, ((childItr-1)*realRFSize+1):childItr*realRFSize, :) = img;
        end
        if realRFSize < 100
              childImgs = imresize(childImgs, 2);
        end
        imwrite(childImgs, [debugFolder '/part' num2str(vocabItr) '.png']);
        
        %% Assign distributions.
        nodeDistributionLevel(vocabItr).childrenLabelDistributions = childrenLabelDistributions;
        nodeDistributionLevel(vocabItr).childrenPosDistributionProbs = childrenProbs;
        
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