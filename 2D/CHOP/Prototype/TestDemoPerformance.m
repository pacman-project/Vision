function [] = TestDemoPerformance(datasetName, mode)
    % Recognition-related stuff.
    global vocabulary, global redundantVocabulary, global modes, global highLevelModes;
    numberOfTestImages = 250;
    
    %% Read vocabulary.
    % Set program options, and create relevant folders.
    options = SetParameters(datasetName, false);
    createFolders(options);
    
    % Read vocabulary from the file.
    if exist([pwd '/output/' datasetName '/vb.mat'], 'file')
        load([pwd '/output/' datasetName '/vb.mat'], 'vocabulary', 'redundantVocabulary', 'modes', 'highLevelModes');
    else
        display('No vocabulary exists!');
    end
    
    % Read models for predictions from file
    if exist([pwd '/demo/Category_Pose/models/' datasetName '_models.mat'], 'file')
        load([pwd '/demo/Category_Pose/models/' datasetName '_models.mat']);
    else
        display('No prediction models exist!');
        models = [];
    end
    options.prediction.models = models;
    options.prediction.feature_params = feature_params;
    
    %% Main loop.
    % Process input image, and show the results in GUI.
    categoryStrs = {'Mug', 'Bowl', 'Tool', 'Pan'};
    tp = 0;
    fp = 0;
    fn = 0;
    
    % Randomize order of test images to get a glimpse of the performance
    % earlier on.
    if ~exist([datasetName 'randIdx.mat'], 'file')
        % Read all training images, and test them.
        if strcmp(mode, 'classify')
            testFileNames = fuf([pwd '/input/' datasetName '/test/*.png'], 1, 'detail');
        else
            testFileNames = fuf([pwd '/input/' datasetName '/wideTest/*.png'], 1, 'detail');
        end
        randIdx = randperm(numel(testFileNames), numberOfTestImages);
        testFileNames = testFileNames(randIdx);
        save([datasetName 'randIdx.mat'], 'randIdx', 'testFileNames');
    else
        load([datasetName 'randIdx.mat'], 'randIdx', 'testFileNames');
    end

    % Process each image to find detections.
    realCategoryArr = zeros(numel(testFileNames),1);
    detectedArr = zeros(numel(testFileNames),1);
    for fileItr = 222:numel(testFileNames)
        
        img = imread(testFileNames{fileItr});
        partImg = processImage(img, datasetName, testFileNames{fileItr}, options);
        label = setdiff(unique(partImg),0);
        
%         % Get real category of the image.
%         for categoryItr = 1:numel(categoryStrs)
%             isCorrectCategory = strfind(testFileNames{fileItr}, categoryStrs{categoryItr});
%             if ~isempty(isCorrectCategory)
%                 realCategoryArr(fileItr) = categoryItr;
%                 break;
%             end
%         end
%         
%         if isempty(label)
%             fn = fn + 1;
%             detectedArr(fileItr) = -1;
%             continue;
%         end
%             
%         detections = setdiff(unique(partImg),0);
%         detected = false;
%         for detItr = 1:numel(detections)
%             if strfind(testFileNames{fileItr}, categoryStrs{detections(detItr)})
%                 % That means, the category has been found right.
%                 tp = tp + 1;
%                 detected = true;
%                 detectedArr(fileItr) = detections(detItr);
%             else
%                 fp = fp + 1;
%                 if ~detected
%                     detectedArr(fileItr) = detections(detItr);
%                 end
%             end
%         end
%         
%         % Get real category of the image.
%         for categoryItr = 1:numel(categoryStrs)
%             isCorrectCategory = strfind(testFileNames{fileItr}, categoryStrs{categoryItr});
%             if ~isempty(isCorrectCategory)
%                 realCategoryArr(fileItr) = categoryItr;
%                 break;
%             end
%         end
%         
%         % If not found, increase false negative count.
%         if ~detected
%             fn = fn +1;
%         end
%         
%         precision = tp / (tp + fp);
%         recall = tp / (tp + fn);
%         display([datasetName ' performance, precision: ' num2str(precision) ', recall: ' num2str(recall) '.']);
%         
%         % Save accuracy array from time to time.
%         if mod(fileItr, 10) == 0 || fileItr == numel(testFileNames)
%             save([datasetName '_acc.mat'], 'realCategoryArr', 'detectedArr');
%         end
    end
    save([pwd '/' datasetName '_perf.mat'], 'tp', 'fp', 'fn');
    precision = tp / (tp + fp);
    recall = tp / (tp + fn);
    
    display([' Final ' datasetName ' performance, precision: ' num2str(precision) ', recall: ' num2str(recall) '.']);
end

% %% Process the image and get level 4 part detections.
% function partImg = processImage(img, datasetName, testFileName, options)
%     % Copy file to the input folder for test.
%     inputImageName = 'lastFrameTestDemo';
%     if ~exist([pwd '/input/' datasetName '/test'], 'dir')
%         mkdir([pwd '/input/' datasetName '/test']);
%     end
%     imwrite(img, [pwd '/input/' datasetName '/test/' inputImageName '.png']);
%     
%     % Run test, and return the output image, along with realizations.
%     [~, fileName, ~] = fileparts(testFileName);
%     if ~exist([pwd '/output/' datasetName '/test/inference/' fileName '_test.mat'], 'file')
%         singleTestImage(testFileName, options);
%     end
%     display(['File name:' fileName]);
%     img = imread([pwd '/output/' datasetName '/smoothed/' fileName '.png']);
%     load([pwd '/output/' datasetName '/test/inference/' fileName '_test.mat'], 'exportArr');
%     
%     % Update feature parameters.
%     options.prediction.feature_params.isTesting = 1;
%     
%     % Detect objects in the scene.
%     [detections, windowSizes] = FindDetections(exportArr, options.prediction.models, options.prediction.feature_params, [size(img,1), size(img,2)]);
%     overlayImg = zeros(size(img,1), size(img,2));
%     if isempty(detections)
%         partImg = overlayImg;
%         return;
%     end
%     
%     % Visualize detection frames.
%     for detItr = 1:size(detections,1)
%         halfSizes = floor((windowSizes{detections(detItr, 3)} / 2) - 3);
%         overlayImg((detections(detItr, 1)-halfSizes(1)):(detections(detItr, 1)+(halfSizes(1)-1)), ...
%             [(detections(detItr, 2)-halfSizes(2)), (detections(detItr, 2)+(halfSizes(2)-1))]) = detections(detItr,4);
%         overlayImg([(detections(detItr, 1)-halfSizes(1)),(detections(detItr, 1)+(halfSizes(1)-1))], ...
%             (detections(detItr, 2)-halfSizes(2)):(detections(detItr, 2)+(halfSizes(2)-1))) = detections(detItr,4);
%     end
%     
%     partImg = overlayImg;
%     img = imresize(img, size(partImg));
%     img(partImg>0) = 255;
%     imwrite(img, [pwd '/demo/output/' fileName '.png']);
% end
%% Process the image and get level 4 part detections.
function partImg = processImage(img, datasetName, testFileName, options)
    % Copy file to the input folder for test.
    categoryStrs = {'Mug', 'Bowl', 'Tool', 'Pan'};
    inputImageName = 'lastFrameTestDemo';
    if ~exist([pwd '/input/' datasetName '/test'], 'dir')
        mkdir([pwd '/input/' datasetName '/test']);
    end
    imwrite(img, [pwd '/input/' datasetName '/test/' inputImageName '.png']);
    
    [~, fileName, ~] = fileparts(testFileName);
    delete([pwd '/output/' datasetName '/test/inference/' fileName '_test.mat']);
    
    % Run test, and return the output image, along with realizations.
    singleTestImage(testFileName, options);
    img = imread([pwd '/output/' datasetName '/smoothed/' fileName '.png']);
    orgImg = imread([pwd '/output/' datasetName '/original/' fileName '.png']);
    
    if exist([pwd '/output/' datasetName '/test/inference/' fileName '_test.mat'], 'file')
      load([pwd '/output/' datasetName '/test/inference/' fileName '_test.mat'], 'exportArr');
    else
      exportArr = [];
    end
    
    % Update feature parameters.
    options.prediction.feature_params.isTesting = 1;
    overlayImg = zeros(size(img,1), size(img,2));
    
    % If no nodes exist, exit.
    if isempty(exportArr)
        partImg = overlayImg;
        return;
    end
    
    % Detect objects in the scene.
    [detections, ~, groupMask] = FindDetections(exportArr, options.prediction.models, options.prediction.feature_params, [size(img,1), size(img,2)]);
    
    if isempty(detections)
        partImg = overlayImg;
        return;
    end
    
    % Visualize detection frames.
    detMask = zeros(size(img,1), size(img,2));
    objectVisImg = zeros(size(orgImg), 'uint8');
    for detItr = 1:size(detections,1)
        detMask(detections(detItr, 1), detections(detItr, 2)) = detections(detItr,4);
        
        % Visualize results.
        objectMaskStr = [options.currentFolder '/render/' categoryStrs{detections(detItr,4)} ...
            '/' categoryStrs{detections(detItr,4)} num2str(round(detections(detItr,5)/30)*30) '.png'];
        objectMask = imread(objectMaskStr);
        binaryObjectMask = rgb2gray(objectMask)>0;
        objectMaskSize = [size(objectMask,1), size(objectMask,2)];
        lowBounds = detections(detItr, 1:2) - round([size(objectMask,1), size(objectMask,2)]/2);
        for bandItr = 1:size(orgImg,3)
            overlapImg = objectVisImg(lowBounds(1):(lowBounds(1)+objectMaskSize(1)-1), lowBounds(2):(lowBounds(2)+objectMaskSize(2)-1), bandItr);
            objectMaskBand = objectMask(:,:,bandItr);
            overlapImg(binaryObjectMask) = objectMaskBand(binaryObjectMask);
            objectVisImg(lowBounds(1):(lowBounds(1)+objectMaskSize(1)-1), lowBounds(2):(lowBounds(2)+objectMaskSize(2)-1), bandItr) = overlapImg;
        end
    end
    
    % TODO: Show objectVisImg in gui.
    partImg = objectVisImg;
    
    % Get highest-level clean reconstructed image and multiply it with
    % group image. Each object's realizations will thus have a different
    % colour.
    maxLevelId = min(1, (max(exportArr(:,4))-2));
    cleanRealizationMask = double(imread([options.outputFolder '/reconstruction/test/' fileName '_level' num2str(maxLevelId) 'clean.png']));
    cleanRealizationMask = cleanRealizationMask / max(max(cleanRealizationMask));
    groupImg = double(label2rgb(bwlabel(groupMask), 'jet', 'k', 'shuffle'));
    for bandItr = 1:size(groupImg,3)
        groupImg(:,:,bandItr) = groupImg(:,:,bandItr) .* cleanRealizationMask;
    end
    groupImg = uint8(round(groupImg));
    
    % Combine group image with smoothed image.
    cleanRealizationMask = uint8(cleanRealizationMask>0);
    for bandItr = 1:size(groupImg,3)
        groupImgBand = groupImg(:,:,bandItr);
        groupImg(:,:,bandItr) = groupImgBand .* cleanRealizationMask + img .* uint8(~cleanRealizationMask);
    end
    
    % TODO: Show groupImg in gui.
    
%     
%     % Visualize detection frames.
%     for detItr = 1:size(detections,1)
%         halfSizes = floor((windowSizes{detections(detItr, 3)} / 2) - 3);
%         overlayImg((detections(detItr, 1)-halfSizes(1)):(detections(detItr, 1)+(halfSizes(1)-1)), ...
%             [(detections(detItr, 2)-halfSizes(2)), (detections(detItr, 2)+(halfSizes(2)-1))]) = detections(detItr,4);
%         overlayImg([(detections(detItr, 1)-halfSizes(1)),(detections(detItr, 1)+(halfSizes(1)-1))], ...
%             (detections(detItr, 2)-halfSizes(2)):(detections(detItr, 2)+(halfSizes(2)-1))) = detections(detItr,4);
%     end
%     
%     partImg = overlayImg;
%     img = imresize(img, size(partImg));
%     img(partImg>0) = 255;
    imwrite(orgImg, [pwd '/demo/output/' fileName '_0_org.jpg']);
    imwrite(img, [pwd '/demo/output/' fileName '_1_smoothed.png']);
    imwrite(partImg, [pwd '/demo/output/' fileName '_2_detection.png']);
    imwrite(groupImg, [pwd '/demo/output/' fileName '_3_parts.png']);
end