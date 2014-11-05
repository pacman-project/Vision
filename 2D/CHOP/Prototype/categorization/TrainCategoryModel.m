function [] = TrainCategoryModel(datasetName, minLevel, maxLevel, poolSize, imgSize)
    if ~exist([pwd '/models/' datasetName '_data_' num2str(minLevel) '_' num2str(maxLevel) '.mat'], 'file')
        % Load relevant info.
        load([pwd '/output/' datasetName '/vb.mat']);
        load([pwd '/output/' datasetName '/export.mat']);
    
        if maxLevel > numel(vocabulary)
            maxLevel = numel(vocabulary);
        end
        if minLevel  < 1 || minLevel>maxLevel
            minLevel = 1;
        end
        featureDims = cellfun(@(x) numel(x), vocabulary);
        % Extract simple features from images.
        cumSums = int32(cumsum(featureDims));
        cumSums(2:end) = cumSums(1:(end-1));
        cumSums(1) = 0;
        cumSums(end+1) = sum(featureDims);
        cumSums = cumSums(:);
        featureDim = sum(featureDims(minLevel:maxLevel));
        numberOfImages = max(exportArr(:,5));
        features = zeros(numberOfImages, featureDim * poolSize * poolSize);
        allActivations = exportArr;
        allActivations(:,1) = allActivations(:,1) + cumSums(allActivations(:,4));
        allActivations = allActivations(allActivations(:,4) >= minLevel & allActivations(:,4) <= maxLevel,:);
        allActivations(:,1) = allActivations(:,1)- cumSums(minLevel);
%        allActivations = unique(allActivations(:,[1 3]), 'rows');
        poolSteps = imgSize / poolSize;
        for imgItr = 1:numberOfImages
            activations = allActivations(allActivations(:,5) == imgItr,:);
            for poolItr1 = 1:poolSize
                for poolItr2 = 1:poolSize
                    minX = round(((poolItr1-1) * poolSteps(1)+1));
                    maxX = round(((poolItr1) * poolSteps(1)));
                    minY = round(((poolItr2-1) * poolSteps(2)+1));
                    maxY = round(((poolItr2) * poolSteps(2)));
                    newActivations  = activations(activations(:,2) >= minX & activations(:,2) <= maxX & activations(:,3) >= minY & activations(:,3) <= maxY,:);
                    startIdx = ((poolItr1-1) * poolSize + (poolItr2-1))*featureDim + 1;
                    tempFeatureVect = features(imgItr,startIdx:(startIdx+featureDim-1));
                    for actItr = 1:size(newActivations,1)
                        tempFeatureVect(newActivations(actItr,1)) = tempFeatureVect(newActivations(actItr,1)) + 1;
                    end
                    features(imgItr,startIdx:(startIdx+featureDim-1)) = tempFeatureVect;
                end
            end
        end
        normFactor = max(max(features));
        features = features / max(max(features));
        features = double(features > 0);
        W = [];
        
        %% Apply linear discriminant analysis
        % First, prevent overfitting by applying PCA.
%         [COEFF, scores] = princomp(features);
%         scores = scores(:,1:20);
%         [features, W] = FDA(scores', categoryArrIdx);
%         features = features';
        % Train a SVM model with cross-validation.
        %% parameter selection for 1-vs-1 multi-class classification
         bestcv = 0;
         bestc = -1;
         bestg = -1;
   %      for log2c = 7:15,
   %        for log2g = [1/32, 1/16, 1/8,1/4,  1/2, 1, 2, 4, 8]
        for log2c = 7
            for log2g = 1/32;
             cmd = ['-v 10 -t 2 -c ', num2str(log2c), ' -g ', num2str(log2g),' -q '];
             cv = svmtrain(categoryArrIdx, features, cmd);
             if (cv >= bestcv),
               bestcv = cv; bestc = log2c; bestg = log2g;
             end
             fprintf('%g %g %g (best c=%g, g=%g, rate=%g)\n', log2c, log2g, cv, bestc, bestg, bestcv);
           end
         end
        cmd = ['-t 2 -c ', num2str(bestc), ' -g ', num2str(bestg),' -q '];
        learnedModel = svmtrain(categoryArrIdx, features, cmd);
        save([pwd '/models/' datasetName '_data_' num2str(minLevel) '_' num2str(maxLevel) '.mat'], 'features', 'categoryArrIdx', 'learnedModel', 'cumSums', 'W', 'normFactor', 'COEFF');
    end
end