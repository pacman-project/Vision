function [validSubs, overallCoverage, overallMatchCost] = getMRMRParts(bestSubs, numberOfFinalSubs, ...
    categoryArrIdx, imageIdx, validationIdx, valItr, midThr, singlePrecision)

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
    
    % Learn train/validation features, and assign labels.
    for bestSubItr = 1:numberOfBestSubs
        instanceImageIdx = allInstances(allInstances(:,2) == bestSubItr, 1);
        
        % For every instance of every sub in every image, we add 1 to the 
        % corresponding counter location.
        for instanceItr = 1:numel(instanceImageIdx)
            allFeatures(instanceImageIdx(instanceItr), bestSubItr) = ...
                allFeatures(instanceImageIdx(instanceItr), bestSubItr) + 1;
        end
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
    validSubs = sort(validSubs);
    
    % Get train and validation accuracy for final evaluation.
    [trainAccuracy, ~] = calculateCategorizationAccuracy(bestSubs(validSubs), ...
       categoryArrIdx, imageIdx, validationIdx, valItr, midThr, singlePrecision, 0, true);
    [trueAccuracy, ~] = calculateCategorizationAccuracy(bestSubs(validSubs), ...
       categoryArrIdx, imageIdx, validationIdx, valItr, midThr, singlePrecision, 1, true);
    display(['[SUBDUE] [Val ' num2str(valItr) '] Val accuracy after discriminative part selection : %' ...
        num2str(100 * trueAccuracy) ...
        ', with training set accuracy : %' num2str(100 * trainAccuracy), ...
        ', having ' num2str(numel(validSubs)) ' subs.']);  
end