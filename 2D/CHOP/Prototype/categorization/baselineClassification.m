function [ ] = baselineClassification(datasetName)
    %BASELINEFEATUREEXTRACTION Summary of this function goes here
    poolSizes = [1, 2, 3];
    
    %   Detailed explanation goes here
    options = SetParameters(datasetName, true);
    datasetTestFolder = [options.currentFolder '/output/' datasetName '/test/inference/'];
    load([options.currentFolder '/output/' datasetName '/vb.mat'], 'vocabulary', 'categoryNames');
    load([options.currentFolder '/output/' datasetName '/export.mat'], 'exportArr', 'categoryArrIdx', 'trainingFileNames');
    testFileNames = fuf([datasetTestFolder '*.mat'], 1, 'detail');
    
    %% Step 1.1: Extract a set of features from the input images.
    display('..... Level 1 Feature Extraction started. This may take a while.');
    allFeatures = cell(numel(trainingFileNames), numel(vocabulary), numel(poolSizes));
    for fileItr = 1:numel(trainingFileNames)
        
        % Get the size of the image.
        img = imread(trainingFileNames{fileItr});
        imgSize = size(img);
        if numel(imgSize) > 2
            imgSize = imgSize(1:2);
        end
        imgArr = exportArr(exportArr(:,5) == fileItr, :);
        
        for levelItr = 1:numel(vocabulary)
            vocabLevel = vocabulary{levelItr};
            levelArr = imgArr(imgArr(:,4) == levelItr,:);
            for poolSizeItr = poolSizes
                poolSize = poolSizes(poolSizeItr);
                stepSizes = ceil(imgSize/poolSize);
                
                imgFeatures = zeros(1, numel(vocabLevel) * poolSize * poolSize);
                
                for poolItr1 = 1:poolSize
                    for poolItr2 = 1:poolSize

                        minX = (poolItr1-1) * stepSizes(1);
                        minY = (poolItr2-1) * stepSizes(2);
                        maxX = min((imgSize(1)+1), (poolItr1 * stepSizes(1)));
                        maxY = min((imgSize(2)+1), (poolItr2 * stepSizes(2)));

                        % Get the features that belong to this bin.
                        poolArr = levelArr(levelArr(:,2) >= minX & levelArr(:,2) < maxX & ...
                                    levelArr(:,3) >= minY & levelArr(:,3) < maxY, 1);
                        startOffset = (poolItr1-1) * poolSize * numel(vocabLevel) + (poolItr2-1) * numel(vocabLevel);
                        for featureItr = 1:size(poolArr,1)
                            imgFeatures(startOffset + poolArr(featureItr)) = imgFeatures(poolArr(featureItr)) + 1;
                        end
                    end
                end
                
                % Normalize feature vector.
                if nnz(imgFeatures) > 0
                    imgFeatures = normr(imgFeatures);
                end
                allFeatures{fileItr, levelItr, poolSizeItr} = imgFeatures;
            end
        end
    end
    
    %% Here, we will run training using the features learned.
    % A separate model is learned for every pooling & vocabulary level pair.
    
    for levelItr = 1:numel(vocabulary)
        for poolSizeItr = 1:numel(poolSizes)
            relevantFeatures = cat(1, allFeatures{:, levelItr, poolSizeItr});
            validRows = sum(relevantFeatures,2) > 0;
            trainLabels = categoryArrIdx(validRows, :);
            trainFeatures = relevantFeatures(validRows, :);
            
            bestcv = 0;
            bestc = -1;
            for log2c = [1/64, 1/32, 1/16, 1/8, 1/4, 1/2, 1, 2, 4, 8, 16, 32, 64, 128]
                cmd = ['-v 5 -t 0 -c ', num2str(log2c),' -q '];
                cv = svmtrain(trainLabels, trainFeatures, cmd);
                if (cv >= bestcv),
                    bestcv = cv; bestc = log2c;
                end
                fprintf('%g %g (best c=%g, rate=%g)\n', log2c, cv, bestc, bestcv);
            end
            cmd = ['-t 0 -c ', num2str(bestc), ' -q '];
            learnedModel = svmtrain(trainLabels, trainFeatures, cmd);
            save([pwd '/models/' datasetName '_level' num2str(levelItr) '_pool' num2str(poolSizes(poolSizeItr)) '.mat'], 'learnedModel');
        end
    end
    
    %% Collect features for test images.
    testFeatures = cell(numel(testFileNames), numel(vocabulary), numel(poolSizes));
    testLabels = zeros(numel(testFileNames), 1);
    for fileItr = 1:numel(testFileNames)
        
        % Get the size of the image.
        load(testFileNames{fileItr}, 'imgSize', 'exportArr', 'categoryLabel');
        testLabels(fileItr) = categoryLabel;
        
        for levelItr = 1:numel(vocabulary)
            vocabLevel = vocabulary{levelItr};
            levelArr = exportArr(exportArr(:,4) == levelItr,:);
            for poolSizeItr = poolSizes
                poolSize = poolSizes(poolSizeItr);
                stepSizes = ceil(imgSize/poolSize);
                
                imgFeatures = zeros(1, numel(vocabLevel) * poolSize * poolSize);
                
                for poolItr1 = 1:poolSize
                    for poolItr2 = 1:poolSize

                        minX = (poolItr1-1) * stepSizes(1);
                        minY = (poolItr2-1) * stepSizes(2);
                        maxX = min((imgSize(1)+1), (poolItr1 * stepSizes(1)));
                        maxY = min((imgSize(2)+1), (poolItr2 * stepSizes(2)));

                        % Get the features that belong to this bin.
                        poolArr = levelArr(levelArr(:,2) >= minX & levelArr(:,2) < maxX & ...
                                    levelArr(:,3) >= minY & levelArr(:,3) < maxY, 1);
                        startOffset = (poolItr1-1) * poolSize * numel(vocabLevel) + (poolItr2-1) * numel(vocabLevel);
                        for featureItr = 1:size(poolArr,1)
                            imgFeatures(startOffset + poolArr(featureItr)) = imgFeatures(poolArr(featureItr)) + 1;
                        end
                    end
                end
                
                % Normalize feature vector.
                if nnz(imgFeatures) > 0
                    imgFeatures = normr(imgFeatures);
                end
                testFeatures{fileItr, levelItr, poolSizeItr} = imgFeatures;
            end
        end
    end
    
    %% Finally, make predictions and estimate the accuracy.
    for poolSizeItr = 1:numel(poolSizes)
        existingPredLabels = -1 * ones(numel(testFileNames),1);
        for levelItr = 1:numel(vocabulary)
            curFeatures = cat(1, testFeatures{:, levelItr, poolSizeItr});
            validRows = sum(curFeatures,2) > 0;
            load([pwd '/models/' datasetName '_level' num2str(levelItr) '_pool' num2str(poolSizes(poolSizeItr)) '.mat'], 'learnedModel');
            cmd = '-q';
            [predLabels,~, ~] = svmpredict(testLabels, curFeatures, learnedModel, cmd);
            predLabels(~validRows) = -1;
            
            % Combine with existing labels.
            existingPredLabels(predLabels ~= -1) = predLabels(predLabels ~= -1);
            predLabels(predLabels == -1) = existingPredLabels(predLabels==-1);
            accuracy = nnz(predLabels == testLabels) / numel(testLabels); %#ok<NASGU>
            confMat = confusionmat(testLabels, predLabels) %#ok<NASGU,NOPRT>
            % Mark invalid rows. 
            save([pwd '/models/' datasetName '_level' num2str(levelItr) '_pool' num2str(poolSizes(poolSizeItr)) '.mat'], 'accuracy', 'confMat', '-append');
        end
    end
end

function [labels] = knnclassify(trainFeatures, trainLabels, testFeatures, k, distance)
    [IDX, ~] = knnsearch(trainFeatures, testFeatures, 'K', k, 'Distance', distance);
    IDX = trainLabels(IDX);
    featureDim = size(IDX,2);
    for itr = 1:size(IDX,1)
       ordering = randperm(featureDim);
       IDX(itr,:) = IDX(itr, ordering);
       [vals, ~, IC] = unique(IDX(itr,:), 'stable');
       IDX(itr,1) = vals(mode(IC));
    end
    labels = IDX(:,1);
end