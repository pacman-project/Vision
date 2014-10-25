function [ accuracy ] = TestCategoryModel(datasetName, minLevel, maxLevel)
    % Load relevant info.
    load([pwd '/models/' datasetName '_data_' num2str(minLevel) '_' num2str(maxLevel) '.mat']);
%    load([pwd '/models/' datasetName '_testLabels.mat']);
    testData = fuf([pwd '/output/' datasetName '/test/inference/*.mat'], 1, 'detail');
    if numel(cumSums) < maxLevel+1
        maxLevel = numel(cumSums)-1;
    end
    features = zeros(numel(testData), size(learnedModel.SVs,2));
    cumSums = cumSums - cumSums(minLevel);
    testLabels = zeros(numel(testData),1);
    for imgItr = 1:numel(testData);
        load(testData{imgItr});
        testLabels(imgItr) = categoryLabel;
        testActivations = exportArr(:, [1 4 5]);
        testActivations = testActivations(ismember(testActivations(:,2), minLevel:maxLevel),:);
        testActivations(:,1) = testActivations(:,1) + cumSums(testActivations(:,2));
        testActivations = testActivations(:,1);
        features(imgItr, testActivations) = 1;
    end
    features = features(:,(cumSums(minLevel)+1):(cumSums(maxLevel+1)));
    
    cmd='';
    %  [predicted_category_label, accuracy, prob_estimates] = svmpredict(labels_sub, features, LearnedModels.category_model, cmd);
    [predLabels,~, ~] = svmpredict(testLabels, features, learnedModel, cmd);
    accuracy = numel(find(testLabels==predLabels))/ numel(predLabels);
end