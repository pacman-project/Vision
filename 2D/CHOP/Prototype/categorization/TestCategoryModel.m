function [ accuracy ] = TestCategoryModel(datasetName, minLevel, maxLevel)
    % Load relevant info.
    load([pwd '/models/' datasetName '_data_' num2str(minLevel) '_' num2str(maxLevel) '.mat']);
    load([pwd '/models/' datasetName '_testLabels.mat']);
    testData = fuf([pwd '/output/' datasetName '/test/inference/*.mat'], 1, 'detail');
    features = zeros(numel(testData), size(learnedModel.SVs,2));
    for imgItr = 1:numel(testData);
        load(testData{imgItr});
        testActivations = exportArr(:, [1 4 5]);
        testActivations = testActivations(ismember(testActivations(:,2), minLevel:maxLevel),:);
        testActivations(:,1) = testActivations(:,1) + cumSums(testActivations(:,2));
        testActivations = testActivations(:,1);
        features(imgItr, testActivations) = 1;
    end
    features = features(:,(cumSums(minLevel)+1):end);
    if maxLevel ~= numel(cumSums)
       features = features(:,1:(cumSums(maxLevel+1))); 
    end
    
    cmd='';
    %  [predicted_category_label, accuracy, prob_estimates] = svmpredict(labels_sub, features, LearnedModels.category_model, cmd);
    [predLabels,~, ~] = svmpredict(testLabels, features, learnedModel, cmd);
    accuracy = numel(find(testLabels==predLabels))/ numel(predLabels);
end