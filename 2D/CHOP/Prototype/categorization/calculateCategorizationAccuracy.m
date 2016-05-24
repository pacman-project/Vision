function [bestAcc, bestPrecision] = calculateCategorizationAccuracy(bestSubs, ...
    categoryArrIdx, imageIdx, validationIdx, valItr, midThr, singlePrecision, includeValidationSet, optimizeSVM)
    % Initialize data structures.
    categoryArrIdx = categoryArrIdx';
    numberOfBestSubs = numel(bestSubs);
    
    % Put the instances ([labelId, imageId, validationIdx] triple) into an array.
    allInstances = cell(numel(bestSubs),1);
    for bestSubItr = 1:numel(bestSubs)
       adaptiveThreshold = ((midThr * (size(bestSubs(bestSubItr).edges,1) * 2 + 1)) + singlePrecision);   
       allInstances{bestSubItr} = [imageIdx(bestSubs(bestSubItr).instanceCenterIdx),...
           repmat(bestSubItr, numel(bestSubs(bestSubItr).instanceCenterIdx), 1), ...
           bestSubs(bestSubItr).instanceValidationIdx];
    end
    allInstances = cat(1, allInstances{:});
    allInstances = sortrows(allInstances);
    
    % Get image ids for training and validation set.
    trainingImageIdx = unique(allInstances(allInstances(:,3) ~= valItr,1));
    validationImageIdx = unique(imageIdx(validationIdx == valItr));
    
    % Obtain train/validation data labels.
    [uniqueImageIdx, IA, ~] = unique(imageIdx);
    imageCategoryIdx = zeros(max(uniqueImageIdx),1);
    imageCategoryIdx(uniqueImageIdx) = categoryArrIdx(IA);
    
    % Create the feature vectors to be classified.
    allFeatures = zeros(max(imageIdx), numberOfBestSubs);
    trainLabels = imageCategoryIdx(trainingImageIdx);
    validationLabels = imageCategoryIdx(validationImageIdx);
    imageIds = unique(imageIdx);
    
    % Learn train/validation features, and assign labels.
    for imageItr = imageIds'
        instanceImageIdx = allInstances(allInstances(:,1) == imageItr, 2);
        featureArr = hist(double(instanceImageIdx), 1:numberOfBestSubs);
        allFeatures(imageItr, :) = featureArr;
    end 
    
    % Assign train/test features and normalize them (converting them to
    % unit vectors).
    validRows = sum(allFeatures,2) ~= 0;
 %   allFeatures(validRows,:) = normr(allFeatures(validRows,:));
    allFeatures(validRows,:) = allFeatures(validRows,:) > 0;
    trainFeatures = allFeatures(trainingImageIdx,:);
    validTrainingRows = sum(trainFeatures,2) ~= 0;
    validationFeatures = allFeatures(validationImageIdx,:);
    validValidationRows = sum(validationFeatures,2) ~= 0;
    
    % Get only valid rows for training.
    trainFeatures = trainFeatures(validTrainingRows, :);
    trainLabels = trainLabels(validTrainingRows, :);
    
    % Combine training and validation sets to evaluate the learned model.
    if includeValidationSet
%        validationFeatures = [validationFeatures; trainFeatures];
%        validationLabels = [validationLabels; trainLabels];
%        validValidationRows = [validValidationRows; validTrainingRows];
    else
        validationFeatures = trainFeatures;
        validationLabels = trainLabels;
        validValidationRows = validTrainingRows;
    end
    
    % Finally, we classify the validation data and return the performance.
%    bestc = 1;
    bestAcc = 0;
    bestPrecision = 0;
    if optimizeSVM
 %       log2cArr = [1/128, 1/64, 1/32, 1/16, 1/8,1/4, 1/2, 1, 2, 4, 8, 16, 32, 64, 128];
        log2cArr = [1/128, 1/32, 1/8, 1/2, 1, 4, 16, 64, 128];
    else
        log2cArr = 1;
    end
    
    for log2c = log2cArr
        cmd = ['-t 0 -c ', num2str(log2c), ' -q '];
        learnedModel = svmtrain(double(trainLabels), trainFeatures, cmd);
        cmd = '-q';
        [predLabels,~, ~] = svmpredict(double(validationLabels), validationFeatures, learnedModel, cmd);
        predLabels(~validValidationRows) = -1;
        accuracy = nnz(predLabels == validationLabels) / numel(validationLabels);
        precision = nnz(predLabels == validationLabels) / nnz(predLabels~=-1);
        if accuracy > bestPrecision
%            bestc = log2c;
            bestAcc = accuracy;
            bestPrecision = precision;
        end
    end
    
    % Build a class contribution array.
    classContrArr = zeros(numel(bestSubs), numel(unique(trainLabels)));
    for bestSubItr =  1:numel(bestSubs)
        features = trainFeatures(:, bestSubItr) > 0;
        classContrArr(bestSubItr,:) = normr(hist(trainLabels(features), 1:size(classContrArr,2)));
    end
    validationFeatures = validationFeatures>0;
    predLabels = zeros(size(validationLabels));
    for imgItr = 1:numel(predLabels)
       [~, predLabels(imgItr)] = max(sum(classContrArr(validationFeatures(imgItr,:)', :),1));
    end
    bestAcc = nnz(predLabels == validationLabels) / numel(validationLabels);
    bestPrecision = nnz(predLabels == validationLabels) / nnz(predLabels~=-1);
end