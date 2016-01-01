function [ accuracy, confMat ] = TestCategoryModel(datasetName, minLevel, maxLevel, poolSize, imgSize)
    % Load relevant info.
    load([pwd '/models/' datasetName '_data_' num2str(minLevel) '_' num2str(maxLevel) '_' num2str(poolSize) '.mat']);
    load([pwd '/output/' datasetName '/vb.mat'], 'categoryNames');
%    load([pwd '/models/' datasetName '_testLabels.mat']);
    testData = fuf([pwd '/output/' datasetName '/test/inference/*.mat'], 1, 'detail');
    if numel(cumSums) < maxLevel+1
        maxLevel = numel(cumSums)-1;
    end
    poolSteps = imgSize / poolSize;
    featureDim = cumSums(maxLevel+1) - cumSums(minLevel);
    features = zeros(numel(testData), featureDim * (poolSize^2));
    cumSums = cumSums - cumSums(minLevel);
    testLabels = zeros(numel(testData),1);
    for imgItr = 1:numel(testData)
        load(testData{imgItr});
        testLabels(imgItr) = categoryLabel;
        testActivations = exportArr;
        testActivations = testActivations(ismember(testActivations(:,4), minLevel:maxLevel),:);
        testActivations(:,1) = testActivations(:,1) + cumSums(testActivations(:,4));
        for poolItr1 = 1:poolSize
            for poolItr2 = 1:poolSize
                minX = round(((poolItr1-1) * poolSteps(1)+1));
                maxX = round(((poolItr1) * poolSteps(1)));
                minY = round(((poolItr2-1) * poolSteps(2)+1));
                maxY = round(((poolItr2) * poolSteps(2)));
                newActivations  = testActivations(testActivations(:,2) >= minX & testActivations(:,2) <= maxX & testActivations(:,3) >= minY & testActivations(:,3) <= maxY,:);
                startIdx = ((poolItr1-1) * poolSize + (poolItr2-1))*featureDim + 1;
                tempFeatureVect = features(imgItr,startIdx:(startIdx+featureDim-1));
                for actItr = 1:size(newActivations,1)
                    tempFeatureVect(newActivations(actItr,1)) = tempFeatureVect(newActivations(actItr,1)) + 1;
                end
                features(imgItr,startIdx:(startIdx+featureDim-1)) = tempFeatureVect;
            end
        end
    end
    features = normr(features);
%     
    cmd='';
    %  [predicted_category_label, accuracy, prob_estimates] = svmpredict(labels_sub, features, LearnedModels.category_model, cmd);
    [predLabels,~, ~] = svmpredict(testLabels, features, learnedModel, cmd);
    accuracy = numel(find(testLabels==predLabels))/ numel(predLabels);
    
    if ~isempty(strfind(datasetName, 'MNIST'))
        categoryNames = cellfun(@(x) str2double(x), categoryNames);
        customOrder = zeros(size(categoryNames,1),1);
        customOrder(categoryNames+1) = 1:size(categoryNames,1);
    else
        customOrder = [];
    end
    
    confMat = confusionmat(testLabels, predLabels);
    if ~isempty(customOrder)
        confMat= confMat(customOrder, customOrder);
    end
    confMat
end