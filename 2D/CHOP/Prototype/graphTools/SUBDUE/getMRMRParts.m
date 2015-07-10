function [validSubs, overallCoverage, overallMatchCost] = getMRMRParts(bestSubs, numberOfFinalSubs, nodeDistanceMatrix, edgeDistanceMatrix, ...
    categoryArrIdx, imageIdx, validationIdx, valItr, midThr, singlePrecision)

    % Eliminate overlapping subs.
    validSubIdx = find(getDisjointSubs(bestSubs, ...
       nodeDistanceMatrix, edgeDistanceMatrix, singlePrecision, midThr));
    bestSubs = bestSubs(validSubIdx);
    
    if numel(bestSubs) < numberOfFinalSubs
       numberOfFinalSubs = numel(bestSubs); 
    end
    
    % Initialize data structures.
    categoryArrIdx = categoryArrIdx';
    numberOfBestSubs = numel(bestSubs);
    overallCoverage = 0;
    overallMatchCost = 0;
    
    
    % Put the instances ([labelId, imageId, validationIdx] triple) into an array.
    allInstances = cell(numel(bestSubs),1);
    for bestSubItr = 1:numel(bestSubs)
       adaptiveThreshold = ((midThr * (size(bestSubs(bestSubItr).edges,1) * 2 + 1)) + singlePrecision);       
       validInstanceIdx = bestSubs(bestSubItr).instanceMatchCosts < adaptiveThreshold;
       allInstances{bestSubItr} = [imageIdx(bestSubs(bestSubItr).instanceCenterIdx(validInstanceIdx, :)),...
           repmat(bestSubItr, numel(bestSubs(bestSubItr).instanceCenterIdx(validInstanceIdx, :)), 1), ...
           bestSubs(bestSubItr).instanceValidationIdx(validInstanceIdx, :)];
    end
    allInstances = cat(1, allInstances{:});
    allInstances = sortrows(allInstances);
    
    % Get image ids for training and validation set.
    trainingImageIdx = unique(allInstances(allInstances(:,3) ~= valItr,1));
    
    % Obtain train/validation data labels.
    [uniqueImageIdx, IA, ~] = unique(imageIdx);
    imageCategoryIdx = zeros(max(uniqueImageIdx),1);
    imageCategoryIdx(uniqueImageIdx) = categoryArrIdx(IA);
    
    % Create the feature vectors to be classified.
    allFeatures = zeros(max(imageIdx), numberOfBestSubs);
    trainLabels = imageCategoryIdx(trainingImageIdx);
    imageIds = unique(imageIdx);
    
    % Learn train/validation features, and assign labels.
    for imageItr = imageIds'
        instanceImageIdx = allInstances(allInstances(:,1) == imageItr, 2);
        featureArr = hist(instanceImageIdx, 1:numberOfBestSubs);
        allFeatures(imageItr, :) = featureArr;
    end 
    
    % Assign train/test features and normalize them (converting them to
    % unit vectors).
    validRows = sum(allFeatures,2) ~= 0;
    allFeatures(validRows,:) = normr(allFeatures(validRows,:));
    allFeatures = uint8(round(allFeatures * 255));
    trainFeatures = allFeatures(trainingImageIdx,:);
    validTrainingRows = sum(trainFeatures,2) ~= 0;
    
    % Get only valid rows for training.
    trainFeatures = trainFeatures(validTrainingRows, :);
    trainLabels = trainLabels(validTrainingRows, :);
    
    % Here, we apply MR-MR based feature selection.
    validSubs = mrmr_miq_d(trainFeatures, trainLabels, numberOfFinalSubs);
%    validSubs = sort(validSubs);
    
    % Get train and validation accuracy for final evaluation.
    % Finally, we have to select a subset of these features since we would
    % like to maximize the performance with the least number of possible
    % subs.
    maxAcc = 0;
    maxPrec = 0;
    maxAccItr = numel(validSubs);
    if valItr ~= -1
        accArr = zeros(numel(20:20:numel(validSubs)), 1);
        for subItr = 20:20:numel(validSubs)
            [tempAccuracy, tempPrecision] = calculateCategorizationAccuracy(bestSubs(validSubs(1:subItr)), ...
                categoryArrIdx, imageIdx, validationIdx, valItr, midThr, singlePrecision, 1, true);
            accArr(round(subItr/20)) = tempAccuracy;
            if tempPrecision > maxPrec
                maxAccItr = subItr;
                maxPrec = tempPrecision;
%                maxAcc = tempAccuracy;
            end
        end
    end
%    figure, plot(20:20:numel(validSubs), accArr), hold on;
%    axis([20 numel(validSubs) 0 1]);
%    hold off;
    validSubs = validSubs(1:maxAccItr);
    [trueAccuracy, truePrecision] = calculateCategorizationAccuracy(bestSubs(validSubs(1:maxAccItr)), ...
       categoryArrIdx, imageIdx, validationIdx, valItr, midThr, singlePrecision, 1, true);
    display(['[SUBDUE] [Val ' num2str(valItr) '] Val precision after discriminative part selection : %' ...
        num2str(100 * truePrecision) ...
        ', having ' num2str(numel(validSubs)) ' subs.']);  
    
    % Return correct indices now.
    validSubs = validSubIdx(validSubs);
end