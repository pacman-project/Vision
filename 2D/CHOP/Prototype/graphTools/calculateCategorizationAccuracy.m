function [accuracy] = calculateCategorizationAccuracy(bestSubs, categoryArrIdx, imageIdx, validationIdx, valItr)
    % Initialize data structures.
    categoryArrIdx = categoryArrIdx';
    numberOfBestSubs = numel(bestSubs);
    
    % Put the instances ([labelId, imageId, validationIdx] triple) into an array.
    allInstances = cell(numel(bestSubs),1);
    for bestSubItr = 1:numel(bestSubs)
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
    
    % Learn train/validation features, and assign labels.
    for bestSubItr = 1:numberOfBestSubs
        instanceCenterIdx = bestSubs(bestSubItr).instanceCenterIdx;
        instanceImageIdx = imageIdx(instanceCenterIdx);
        
        % For every instance of every sub in every image, we add 1 to the 
        % corresponding counter location.
        for instanceItr = 1:numel(instanceCenterIdx)
            allFeatures(instanceImageIdx(instanceItr), bestSubItr) = ...
                allFeatures(instanceImageIdx(instanceItr), bestSubItr) + 1;
        end
    end 
    
    % Assign train/test features and normalize them (converting them to
    % unit vectors).
    validRows = sum(allFeatures,2) ~= 0;
    allFeatures(validRows,:) = normr(allFeatures(validRows,:));
    trainFeatures = allFeatures(trainingImageIdx,:);
    validTrainingRows = sum(trainFeatures,2) ~= 0;
    validationFeatures = allFeatures(validationImageIdx,:);
    validValidationRows = sum(validationFeatures,2) ~= 0;
    
    % Get only valid rows for training.
    trainFeatures = trainFeatures(validTrainingRows, :);
    trainLabels = trainLabels(validTrainingRows, :);
    
    % Finally, we classify the validation data and return the performance.
    bestc = 1;
    cmd = ['-t 0 -c ', num2str(bestc), ' -q '];
    learnedModel = svmtrain(double(trainLabels), trainFeatures, cmd);
    cmd = '-q';
    [predLabels,~, ~] = svmpredict(double(validationLabels), validationFeatures, learnedModel, cmd);
    predLabels(~validValidationRows) = -1;
    accuracy = nnz(predLabels == validationLabels) / numel(validationLabels);
end