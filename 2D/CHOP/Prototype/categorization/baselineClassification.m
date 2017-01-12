function [ ] = baselineClassification(datasetName)
    %BASELINEFEATUREEXTRACTION Summary of this function goes here
    poolSizes = [1 2];
    
    %   Detailed explanation goes here
    options = SetParameters(datasetName, true);
    datasetTestFolder = [options.currentFolder '/output/' datasetName '/test/inference/'];
    load([options.currentFolder '/output/' datasetName '/vb.mat'], 'vocabulary', 'categoryNames', 'options');
    load([options.currentFolder '/output/' datasetName '/export.mat'], 'categoryArrIdx', 'trainingFileNames');
    if exist([options.currentFolder '/output/' datasetName '/trainingInference.mat'], 'file')
         load([options.currentFolder '/output/' datasetName '/trainingInference.mat'], 'exportArr', 'pooledPositions');
    else
         load([options.currentFolder '/output/' datasetName '/export.mat'], 'exportArr', 'pooledPositions');
         [~, sortIdx] = sort(exportArr(:,5));
         exportArr = exportArr(sortIdx, :);
         pooledPositions = pooledPositions(sortIdx, :);
    end
    testFileNames = fuf([datasetTestFolder '*.mat'], 1, 'detail');
    
    % Set up image dimensions.
    imageSizes = cell(numel(vocabulary), 1);
    if strcmp(options.filterType, 'gabor')
         stride = options.gabor.stride;
    else
         stride = options.auto.stride;
    end
    denom = stride;
    for itr = 1:numel(vocabulary)
         imageSizes(itr) = {ceil(options.imageSize / denom)};
         if ~ismember(itr, options.noPoolingLayers)
              denom = round(denom * options.poolDim);
         end
    end
    
    %% Step 1.1: Extract a set of features from the input images.
    display('..... Level 1 Feature Extraction started. This may take a while.');
    allFeatures = cell(numel(categoryArrIdx), numel(vocabulary), numel(poolSizes));
    maxLevels = numel(vocabulary);
    [~, firstInstances, ~] = unique(exportArr(:, 5), 'R2012a');
    firstInstances = cat(1, firstInstances, size(exportArr,1) + 1);
    exportArr(:, 2:3) = pooledPositions;
    
    for fileItr = 1:numel(categoryArrIdx)
        % Get the size of the image.
        imgArr = exportArr(firstInstances(fileItr):(firstInstances(fileItr+1)-1), :);
      
        for levelItr = 1:numel(vocabulary)
            vocabLevel = vocabulary{levelItr};
            imgSize = imageSizes{levelItr};
            levelArr = imgArr(imgArr(:,4) == levelItr,:);
            for poolSizeItr = 1:numel(poolSizes)
                poolSize = poolSizes(poolSizeItr);
                stepSizes = ceil(imgSize/poolSize);
                
                imgFeatures = zeros(1, numel(vocabLevel) * poolSize * poolSize);
                
                for poolItr1 = 1:poolSize
                    for poolItr2 = 1:poolSize

                        minX = (poolItr1-1) * stepSizes(1) + 1;
                        minY = (poolItr2-1) * stepSizes(2) + 1;
                        maxX = min((imgSize(1)+1), (1 + poolItr1 * stepSizes(1)));
                        maxY = min((imgSize(2)+1), (1 + poolItr2 * stepSizes(2)));

                        % Get the features that belong to this bin.
                        poolArr = levelArr(levelArr(:,2) >= minX & levelArr(:,2) < maxX & ...
                                    levelArr(:,3) >= minY & levelArr(:,3) < maxY, 1);
                        newFeatures = full( sparse( 1, double(poolArr), 1, 1, numel(vocabLevel) ) );
                        startOffset = (poolItr1-1) * poolSize * numel(vocabLevel) + (poolItr2-1) * numel(vocabLevel)+1;
                        imgFeatures(startOffset:(startOffset+(numel(vocabLevel)-1))) = newFeatures;
                    end
                end
                
                % Normalize feature vector.
                if nnz(imgFeatures) > 0
                    imgFeatures = imgFeatures / sum(imgFeatures);
                end
                
                % We apply PCA if needed.
                allFeatures{fileItr, levelItr, poolSizeItr} = imgFeatures;
            end
        end
    end
    
    %% Collect features for test images.
    testFeatures = cell(numel(testFileNames), numel(vocabulary), numel(poolSizes));
    testLabels = zeros(numel(testFileNames), 1);
    for fileItr = 1:numel(testFileNames)
        
        % Get the size of the image.
        load(testFileNames{fileItr}, 'imgSize', 'exportArr', 'categoryLabel', 'pooledPositions');
        exportArr(:, 2:3) = pooledPositions;
        testLabels(fileItr) = categoryLabel;
        
        for levelItr = 1:numel(vocabulary)
            vocabLevel = vocabulary{levelItr};
            imgSize = imageSizes{levelItr};
            levelArr = exportArr(exportArr(:,4) == levelItr,:);
            for poolSizeItr = 1:numel(poolSizes)
                poolSize = poolSizes(poolSizeItr);
                stepSizes = ceil(imgSize/poolSize);
                
                imgFeatures = zeros(1, numel(vocabLevel) * poolSize * poolSize);
                
                for poolItr1 = 1:poolSize
                    for poolItr2 = 1:poolSize

                        minX = (poolItr1-1) * stepSizes(1) + 1;
                        minY = (poolItr2-1) * stepSizes(2) + 1;
                        maxX = min((imgSize(1)+1), (1 + poolItr1 * stepSizes(1)));
                        maxY = min((imgSize(2)+1), (1 + poolItr2 * stepSizes(2)));

                        % Get the features that belong to this bin.
                        poolArr = levelArr(levelArr(:,2) >= minX & levelArr(:,2) < maxX & ...
                                    levelArr(:,3) >= minY & levelArr(:,3) < maxY, 1);
                        newFeatures = full( sparse( 1, double(poolArr), 1, 1, numel(vocabLevel) ) );
                        startOffset = (poolItr1-1) * poolSize * numel(vocabLevel) + (poolItr2-1) * numel(vocabLevel)+1;
                        imgFeatures(startOffset:(startOffset+(numel(vocabLevel)-1))) = newFeatures;
                    end
                end
                
                % Normalize feature vector.
                if nnz(imgFeatures) > 0
                    imgFeatures = imgFeatures / sum(imgFeatures);
                end
                testFeatures{fileItr, levelItr, poolSizeItr} = imgFeatures;
            end
        end
    end
    
    
    %% Here, we will run training using the features learned.
    % A separate model is learned for every pooling & vocabulary level pair.
    for levelItr = 1:maxLevels
        for poolSizeItr = 1:numel(poolSizes)
%              if exist([pwd '/models/' datasetName '_level' num2str(levelItr) '_pool' num2str(poolSizes(poolSizeItr)) '.mat'], 'file')
%                   continue;
%              end
            relevantFeatures = [];
            for itr = 1:levelItr
                 relevantFeatures = cat(2, relevantFeatures, cat(1, allFeatures{:, itr, poolSizeItr}));
            end
 %           relevantFeatures = cat(1, allFeatures{:, levelItr, poolSizeItr});
            validRows = sum(relevantFeatures,2) > 0;
            trainLabels = categoryArrIdx(validRows, :);
            trainFeatures = relevantFeatures(validRows, :);
            
            %% Apply PCA
            colMeans = sum(trainFeatures,1) / size(trainFeatures, 1);
            [coeff, trainFeatures, latent] = princomp(trainFeatures);
            numberOfFeatures = nnz(latent >= 0.002*max(latent));
            trainFeatures = trainFeatures(:, 1:numberOfFeatures);
            
%             %% Grid-search for best parameters.
%             cVals = -1:2:3;
%             gVals = -4:2:0;
%             combs = allcomb(cVals, gVals);
%             cvVals = zeros(size(combs,1),1);
%             parfor itr = 1:size(combs,1)
%                 cmd = ['-v 4 -c ', num2str(2^combs(itr,1)), ' -g ', num2str(2^combs(itr,2)) ' -q'];
%                 cv = svmtrain(trainLabels, trainFeatures, cmd);
%                 fprintf('%g %g %g\n', combs(itr,1), combs(itr,2), cv);
%                 cvVals(itr) = cv;
%             end
%             [~, maxVal] = max(cvVals);
%             cmd = ['-c ', num2str(2^combs(maxVal,1)), ' -g ', num2str(2^combs(maxVal,2)) ' -q'];
            cmd = ['-c ', num2str(2^3), ' -g ', num2str(2^0) ' -q'];
            learnedModel = svmtrain(trainLabels, trainFeatures, cmd);
            save([pwd '/models/' datasetName '_level' num2str(levelItr) '_pool' num2str(poolSizes(poolSizeItr)) '.mat'], 'learnedModel', 'coeff', 'colMeans');
        end
    end
    
    %% Calculate training accuracy.
    for poolSizeItr = 1:numel(poolSizes)
        existingPredLabels = -1 * ones(numel(categoryArrIdx),1);
        for levelItr = 1:maxLevels
             
             
            curFeatures = [];
            for itr = 1:levelItr
                 curFeatures = cat(2, curFeatures, cat(1, allFeatures{:, itr, poolSizeItr}));
            end
%            curFeatures = cat(1, allFeatures{:, levelItr, poolSizeItr});
            validRows = sum(curFeatures,2) > 0;
            load([pwd '/models/' datasetName '_level' num2str(levelItr) '_pool' num2str(poolSizes(poolSizeItr)) '.mat'], 'learnedModel', 'coeff', 'colMeans');
            
            % Apply PCA coefficients.
            curFeatures = curFeatures - repmat(colMeans, size(curFeatures, 1), 1);
            curFeatures = curFeatures * coeff;
            curFeatures = curFeatures(:, 1:size(learnedModel.SVs,2));
            
            cmd = '-q';
            [predLabels,~, ~] = svmpredict(categoryArrIdx, curFeatures, learnedModel, cmd);
            predLabels(~validRows) = -1;
            
            % Combine with existing labels.
            existingPredLabels(predLabels ~= -1) = predLabels(predLabels ~= -1);
            predLabels(predLabels == -1) = existingPredLabels(predLabels==-1);
            
            % Estimate performance, and find confusion matrix.
            trainConfMat = confusionmat(categoryArrIdx, predLabels);
            trainAccuracy = sum(trainConfMat(eye(size(trainConfMat,1)) > 0)) / sum(sum(trainConfMat));
            
            % Mark invalid rows. 
            save([pwd '/models/' datasetName '_level' num2str(levelItr) '_pool' num2str(poolSizes(poolSizeItr)) '.mat'], 'trainAccuracy', 'trainConfMat', '-append');
        end
    end
    
    %% Finally, make predictions and estimate the accuracy.
    for poolSizeItr = 1:numel(poolSizes)
        existingPredLabels = -1 * ones(numel(testFileNames),1);
        for levelItr = 1:maxLevels
             
           curFeatures = [];
            for itr = 1:levelItr
                 curFeatures = cat(2, curFeatures, cat(1, testFeatures{:, itr, poolSizeItr}));
            end
             
 %           curFeatures = cat(1, testFeatures{:, levelItr, poolSizeItr});
            validRows = sum(curFeatures,2) > 0;
            load([pwd '/models/' datasetName '_level' num2str(levelItr) '_pool' num2str(poolSizes(poolSizeItr)) '.mat'], 'learnedModel', 'coeff', 'colMeans');
            cmd = '-q';
            
            % Apply pca coeff.
            curFeatures = curFeatures - repmat(colMeans, size(curFeatures, 1), 1);
            curFeatures = curFeatures * coeff;
            curFeatures = curFeatures(:, 1:size(learnedModel.SVs,2));
            
            [predLabels,~, ~] = svmpredict(testLabels, curFeatures, learnedModel, cmd);
            predLabels(~validRows) = -1;
            
            % Combine with existing labels.
            existingPredLabels(predLabels ~= -1) = predLabels(predLabels ~= -1);
            predLabels(predLabels == -1) = existingPredLabels(predLabels==-1);
            
            % Estimate performance, and find confusion matrix.
            confMat = confusionmat(testLabels, predLabels);
            accuracy = sum(confMat(eye(size(confMat,1)) > 0)) / sum(sum(confMat));
            
            % Mark invalid rows. 
            save([pwd '/models/' datasetName '_level' num2str(levelItr) '_pool' num2str(poolSizes(poolSizeItr)) '.mat'], 'accuracy', 'confMat', '-append');
        end
    end
end

