function [ ] = baselineClassification(datasetName, imageExtension, featureType)
    %BASELINEFEATUREEXTRACTION Summary of this function goes here
    %   Detailed explanation goes here
    poolSize = 3;
    levelItr = 1;
    options = SetParameters(datasetName, true);
    datasetFolder = [options.currentFolder '/input/' datasetName '/vocab/'];
    datasetTestFolder = [options.currentFolder '/input/' datasetName '/test/'];
    load([options.currentFolder '/output/' datasetName '/vb.mat'], 'vocabulary');
    
    %% Step 0.1: Create initial data structures.
    try
        fileNames = fuf([datasetFolder '*', imageExtension], 1, 'detail');
        testFileNames = fuf([datasetTestFolder '*', imageExtension], 1, 'detail');
        fileNames = [fileNames; testFileNames];
        trainTestSep = zeros(numel(fileNames),1)>0;
        trainTestSep((end-(numel(testFileNames)-1)):end) = 1;
    catch
        display(['Unable to find training images. Possibly you forgot to put the images under ./input/' datasetName '/vocab/ or defined the extension wrong.']); 
        return;
    end
    
    %% Learn category/pose of each image and and correct their order 
    % based on the sort order of the images.
    categoryArr = cell(numel(fileNames),1);
    for fileItr = 1:numel(fileNames)
        % Get category label.
        fullName = fileNames{fileItr};
        [~, fileName, ext] = fileparts(fileNames{fileItr});
        if ~isempty(strfind(fileNames{fileItr}, '/vocab/'))
            strLength = numel([options.currentFolder '/input/' options.datasetName '/vocab/']);
        else
            strLength = numel([options.currentFolder '/input/' options.datasetName '/test/']);
        end
        fileNameLength = numel(ext) + numel(fileName) + 1; % 1 for '/' (folder delimiter)
        if numel(fullName) >= strLength + fileNameLength
            categoryStr = fullName(1, (strLength+1):(end-fileNameLength));
            categoryStrSepIdx = strfind(categoryStr, '/');
            if ~isempty(categoryStrSepIdx)
                categoryStr = categoryStr(:, 1:(categoryStrSepIdx(1)-1));
            end
            categoryArr(fileItr) = {categoryStr};
        else
            categoryArr(fileItr) = {''};
        end
    end
    [~, categoryNames, categoryArrIdx] = unique(categoryArr, 'stable');  %#ok<ASGLU>
    categoryNames = categoryArr(categoryNames);
    
    %% Step 1.1: Extract a set of features from the input images.
    if ~exist([pwd '/models/' datasetName '_baseData_.mat'], 'file')
        display('..... Level 1 Feature Extraction started. This may take a while.');
        if strcmp(featureType, 'baseline')
            features = cell(size(fileNames,1),1);
            filters = options.filters;
            filterMatrix = options.filterMatrix;
            if ~exist([options.currentFolder '/categorization/output/' datasetName], 'dir')
               mkdir([options.currentFolder '/categorization/output/' datasetName]);
            end
            for fileItr = 1:size(fileNames,1)
                display(['Processing file ' num2str(fileItr) ' out of ' num2str(numel(fileNames)) ' images.']);
                img = imread(fileNames{fileItr});
                [~, fileName, ~] = fileparts(fileNames{fileItr});

                % Get the Level 1 features.
                [features{fileItr}, responseImg] = getBaselineFeatures(img, filters, filterMatrix, poolSize, options);
                imwrite(responseImg, [options.currentFolder '/categorization/output/' datasetName '/' fileName '.png']);
            end
            features = cat(1, features{:});
        else
            features = getHierarchicalFeatures(datasetName, fileNames, categoryNames, categoryArrIdx, poolSize, numel(vocabulary{levelItr}), levelItr);
        end

        %% Train a svm with these features.
        trainFeatures = features(trainTestSep==0,:);
        trainLabels = categoryArrIdx(trainTestSep==0,:);

        testFeatures = features(trainTestSep==1,:); %#ok<NASGU>
        testLabels = categoryArrIdx(trainTestSep==1,:);
        save([pwd '/models/' datasetName '_baseData_.mat'], 'trainFeatures', 'trainLabels', 'testFeatures', 'testLabels', 'categoryNames');
    else
        load([pwd '/models/' datasetName '_baseData_.mat'], 'trainFeatures', 'trainLabels', 'testFeatures', 'testLabels', 'categoryNames');
    end
        
%     if ~exist([pwd '/models/' datasetName '_baseModel_.mat'], 'file')
%         bestcv = 0;
%         bestc = -1;
%         bestg = -1;
%         for log2c = [1/32, 1/16, 1/8,1/4, 1/2, 1, 2, 4, 8, 16, 32]
%    %     for log2c = -1:6,
%             for log2g = [1/32, 1/16, 1/8,1/4,  1/2, 1, 2, 4, 8]
%       %  for log2c = 7
%      %       for log2g = 1/32;
%                 cmd = ['-v 5 -t 0 -c ', num2str(log2c), ' -g ', num2str(log2g),' -q '];
%                 cv = svmtrain(trainLabels, trainFeatures, cmd);
%                 if (cv >= bestcv),
%                     bestcv = cv; bestc = log2c; bestg = log2g;
%                 end
%                 fprintf('%g %g %g (best c=%g, g=%g, rate=%g)\n', log2c, log2g, cv, bestc, bestg, bestcv);
%             end
%         end
%         cmd = ['-t 0 -c ', num2str(bestc), ' -g ', num2str(bestg),' -q '];
%         learnedModel = svmtrain(trainLabels, trainFeatures, cmd);
%         save([pwd '/models/' datasetName '_baseModel_.mat'], 'learnedModel');
%     else
%         load([pwd '/models/' datasetName '_baseModel_.mat'], 'learnedModel');
%     end
    %% Run classification.
 %   cmd = '-q';
 %   [predLabels,~, ~] = svmpredict(testLabels, testFeatures, learnedModel, cmd);
    predLabels = knnclassify(trainFeatures, trainLabels, testFeatures, 5, 'euclidean');
    accuracy = numel(find(testLabels==predLabels))/ numel(predLabels)
    confMat = confusionmat(testLabels, predLabels) %#ok<NOPRT>
    save([pwd '/models/' datasetName '_baseAccuracy_.mat'], 'accuracy', 'confMat');
end

function [features] = getHierarchicalFeatures(datasetName, fileNames, categoryNames, categoryArrIdx, poolSize, featureCount, levelItr)
    features = cell(size(fileNames,1),1);
    for fileItr = 1:numel(fileNames)
        [~, fileName, ~] = fileparts(fileNames{fileItr});
        img = imread([pwd '/output/' datasetName '/original/' categoryNames{categoryArrIdx(fileItr)} '_' fileName '.png']);
        targetFileName = [pwd '/output/' datasetName '/test/inference/' categoryNames{categoryArrIdx(fileItr)} '_' fileName '_test.mat'];
        load(targetFileName);
        relevantEntries = exportArr(exportArr(:,4) == levelItr, :);
        imgFeatures = zeros(1, featureCount * 4);
        imgHalfSizes = fix(size(img)/poolSize);
        startIdx = 1;
        for poolItr = 1:poolSize
            for poolItr2 = 1:poolSize
               maskedLabels = relevantEntries(relevantEntries(:,2) >= (1+(poolItr-1)*imgHalfSizes(1)) & ...
                   relevantEntries(:,2) < (imgHalfSizes(1)*poolItr) & ...
                   relevantEntries(:,3) >= (1+(poolItr2-1)*imgHalfSizes(2)) & ...
                   relevantEntries(:,3) < (imgHalfSizes(2)*poolItr2), 1);
               imgFeatures(1, startIdx:(startIdx + (featureCount-1))) = hist(maskedLabels, 1:featureCount);
               startIdx = startIdx + featureCount;
            end
        end
        imgFeatures = imgFeatures / sum(imgFeatures);
        features(fileItr) = {imgFeatures};
    end
    features = cat(1, features{:});
end

function [imgFeatures, responseImg] = getBaselineFeatures(img, filters, filterMatrix, poolSize, options)
    
    filterCount = numel(filters);
    responseImgs = zeros(size(img,1), size(img,2), filterCount);
    if strcmp(options.filterType, 'gabor')
        for filtItr = 1:filterCount
            if size(img,3)>1
                img = rgb2gray(img(:,:,1:3));
            end

            %% Apply denoising to get better responses.
            for bandItr = 1:size(img,3);
                img(:,:,bandItr)=medfilt2(img(:,:,bandItr), [3,3]);
            end
            img = double(img);
            currentFilter = double(options.filters{filtItr});
            responseImg = conv2(img, currentFilter, 'same');

            % Save response for future processing
            responseImgs(:,:,filtItr) = responseImg;
        end
    else
        stride = options.auto.stride;
        whMat = options.auto.whMat;
        mu = options.auto.mu;

        filterSize = size(filters{1});
        filterBandSize = filterSize(1:2);
        img = double(img);
        dim1 = (size(img,1)-filterSize(1)+1);
        dim2 = (size(img,2)-filterSize(2)+1);

        % Implementing stride here. Some blocks will be skipped. 
        imgCols = zeros( ceil(dim1/stride) * ceil(dim2/stride), prod(filterSize));
        startIdx = 1;
        iterator = prod(filterBandSize)-1;

        % Get linear indices of correct columns, since not all blocks (each
        % corresponds a column) may be used, depending on the stride parameter.
        idx1 = 1:stride:dim1;
        idx2 = 1:stride:dim2;
        [p,q] = meshgrid(idx1, idx2);
        pairs = [p(:) q(:)];
        pairs = sortrows(pairs,2);
        validCols = sub2ind([dim1, dim2], pairs(:,1), pairs(:,2));

        for bandItr = 1:size(img,3)
            tempCols = im2col(img(:,:,bandItr), filterBandSize)';
            imgCols(:,startIdx:(startIdx+iterator)) = tempCols(validCols,:);
            startIdx = startIdx + iterator + 1;
        end
        clear tempCols;
        muArr = repmat(mu, [size(imgCols,1), 1]);
        halfSize = ceil(filterSize(1)/2);

        % Pre-process the image blocks by whitening them.
        imgCols2 = imgCols - muArr;
        imgCols2 = imgCols2 * whMat;

        % Instead of convolving imgCols2 with the filter, we simply
        % assign a label(filter id) to each row in imgCols2 by
        % finding the cluster with minimum distance to each sample, 
        % given in every row. We also allow a soft-competition of
        % filters by removing roughly half of them, based on the
        % average distance of cluster centers from each sample.
        numberOfCols = size(imgCols2,1);
        % If number of cols is more than 1000, we divide and conquer
        % imgCols2 array into fixed-size bins. Otherwise, it's hell of a
        % memory problem.
        if numberOfCols>1000 % Depends on feature dimension, but 1000 should 
                             % be reasonable for features smaller than 15x15x3)
            numberOfSets = ceil(numberOfCols/1000);
            smallRepFilterMatrix = repmat(filterMatrix, 1000, 1);
            colFiltDistances = zeros(numberOfCols * filterCount,1);
            for setItr = 1:numberOfSets
                colsInSet = min(1000, (numberOfCols-((setItr-1)*1000)));
                setCols = imgCols2((((setItr-1) * 1000) + 1):min(setItr*1000, numberOfCols),:);
                setCols = setCols(floor((0:((colsInSet * filterCount)-1)) / filterCount) + 1, :);
                if colsInSet ~= 1000
                    smallRepFilterMatrix = smallRepFilterMatrix(1:(colsInSet * filterCount),:);
                end

                % Subtract repFilterMatrix from imgCols2.
                totalAssgnCount = colsInSet * filterCount;
                assgnStartIdx = (1000 * filterCount * (setItr-1)) + 1;
                colFiltDistances(assgnStartIdx:...
                    (assgnStartIdx+totalAssgnCount-1)) = ...
                    sqrt(sum((setCols - smallRepFilterMatrix).^2,2));
            end
            clear smallRepFilterMatrix imgCols2 setCols;
        else
            % No need to divide, just find the distances.
            imgCols2 = imgCols2(floor((0:((numberOfCols * filterCount)-1)) / filterCount) + 1, :);
            repFilterMatrix = repmat(filterMatrix, numberOfCols, 1);

            % Subtract repFilterMatrix from imgCols2.
            colFiltDistances = sqrt(sum((imgCols2 - repFilterMatrix).^2,2));
            clear repFilterMatrix imgCols2;
        end
        clear imgCols;
        % Reshape distances into a NxD array, where N is number of columns (blocks),
        % and D is number of filters.
        distancesPerCol = reshape(colFiltDistances, filterCount, numberOfCols).';

        % Find average of every row.
        colMeans = repmat(mean(distancesPerCol,2), 1, filterCount);

        % Suppress distances for every row which is more than its average
        % distance to every filter.
        meanAssgnIdx = distancesPerCol > colMeans;
        distancesPerCol(meanAssgnIdx) = colMeans(meanAssgnIdx);
        responses = (colMeans - distancesPerCol);

        % Assign responses to the actual image.
        realCoordIdx1 = idx1 + halfSize - 1;
        realCoordIdx2 = idx2 + halfSize - 1;
        responseImgs(realCoordIdx1, realCoordIdx2, :) = reshape(responses, [numel(idx1), numel(idx2), filterCount]);
    end
    responseImg = mean(responseImgs,3);
    responseImgs(responseImgs<0) = 0;
    responseImg = (responseImg)/(max(responseImg(:)));

    imgFeatures = zeros(1, filterCount * 4);
    %% Extract features from responseImgs.
    startIdx = 1;
    imgHalfSizes = fix(size(responseImgs)/poolSize);
    for poolItr = 1:poolSize
       for poolItr2 = 1:poolSize
           imgFeatures(1, startIdx:(startIdx + (filterCount-1))) = squeeze(sum(sum(responseImgs((1+(poolItr-1)*imgHalfSizes(1)):(imgHalfSizes(1)*poolItr), ...
                                                                                            (1+(poolItr2-1)*imgHalfSizes(2)):(imgHalfSizes(2)*poolItr2), :),1),2));
       startIdx = startIdx + filterCount;
       end
    end
    imgFeatures = imgFeatures/sum(imgFeatures);
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