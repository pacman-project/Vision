function [] = TrainCategoryModel(datasetName, minLevel, maxLevel)
    if ~exist([pwd '/models/' datasetName '_data_' num2str(minLevel) '_' num2str(maxLevel) '.mat'], 'file')
        % Load relevant info.
        load([pwd '/output/' datasetName '/vb.mat']);
        load([pwd '/output/' datasetName '/export.mat']);
    
        % Extract simple features from images.
        featureDims = cellfun(@(x) numel(x), vocabulary);
        cumSums = cumsum(featureDims);
        cumSums(2:end) = cumSums(1:(end-1));
        cumSums(1) = 0;
        featureDim = sum(featureDims);
        numberOfImages = max(exportArr(:,5));
        features = zeros(numberOfImages, featureDim);
        allActivations = exportArr(:,[1,4,5]);
        allActivations(:,1) = allActivations(:,1) + cumSums(allActivations(:,2));
        allActivations = unique(allActivations(:,[1 3]), 'rows');

        for imgItr = 1:numberOfImages
            activations = allActivations(allActivations(:,2) == imgItr,1);
            features(imgItr,activations) = 1;
        end

        features = features(:,(cumSums(minLevel)+1):end);
        if maxLevel ~= numel(cumSums)
           features = features(:,1:(cumSums(maxLevel+1))); 
        end
        % Train a SVM model with cross-validation.
        cmd='-s 0 -t 0';
        learnedModel = svmtrain(categoryArrIdx, features, cmd);
        save([pwd '/models/' datasetName '_data_' num2str(minLevel) '_' num2str(maxLevel) '.mat'], 'features', 'categoryArrIdx', 'learnedModel', 'cumSums');
    end
end