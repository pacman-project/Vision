%> Name: singleTestImages
%>
%> Description: Process each image, run discovery with vocabulary at each
%> level, and find the class each image belongs to. 
%>
%> @param testFileImages The test image names to work on.
%> @param vocabulary
%> @param redundantVocabulary
%> @param modes
%> @param options Program options.
%>
%> @retval totalInferenceTime The amount of time spent in inference.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 19.12.2013
function [categoryLabel, predictedCategory] = singleTestImage(testFileName, vocabulary, vocabularyDistributions, allModes, modeProbs, categoryName, options)
    %% Get the first level nodes.
    % First, downsample the image if it is too big.
    img = imread(testFileName);
    [~, fileName, ~] = fileparts(testFileName);
    % resize image if necessary.
    if max(size(img)) > options.maxImageDim
       img = imresize(img, options.maxImageDim/max(size(img)), 'bilinear'); 
    end
    imwrite(img, [options.processedFolder '/' categoryName '_' fileName '.png']);
    imgSize = size(img);
    if numel(imgSize) > 2
        imgSize = imgSize(1:2); %#ok<NASGU>
    end

    %% Form the first level nodes.
    if exist([options.testInferenceFolder '/' categoryName '_' fileName '_test.mat'], 'file')
         w = warning ('off','all');
         load([options.testInferenceFolder '/' categoryName '_' fileName '_test.mat'], 'exportArr');
         warning(w);
    end
    
    if ~exist('exportAr', 'var')
         [ nodes, ~, nodeActivations, ~, ~, ~ ] = getNodes( img, [], options );
         if isempty(nodes)
             return;
         end

         % Create output folder for images.
         outputImgFolder = [options.currentFolder '/output/' options.datasetName '/reconstruction/test/' fileName];
         if ~exist(outputImgFolder, 'dir')
              mkdir(outputImgFolder);
         end

         % Save smoothed image.
         % Assign nodes their image ids.
         nodes = int32(cell2mat(nodes));
         [exportArr, activationArr] = inferSubs(fileName, img, vocabulary, vocabularyDistributions, allModes, modeProbs, nodes, nodeActivations, options); %#ok<ASGLU,NASGU>
    else
         load([options.testInferenceFolder '/' categoryName '_' fileName '_test.mat'], 'exportArr', 'activationArr');
    end
    maxLevel = max(exportArr(:,4));
    
%     %% Project stuff from top layer.    % Create data structures required for optimization.
%     if ~exist([pwd '/filters/optimizationFilters.mat'], 'file')
%          [rfSizes, visFilters, optimizedFilters, likelihoodLookupTable] = createOptimizationStructures(options, levelItr, true);
%          save([pwd '/filters/optimizationFilters.mat'], 'rfSizes', 'visFilters', 'optimizedFilters', 'likelihoodLookupTable');
%     else
%          load([pwd '/filters/optimizationFilters.mat'], 'visFilters');
%     end
%     validIdx = exportArr(:,4) == maxLevel;
%     maxLevelNodes = exportArr(validIdx,:);
%     maxLevelActivations = activationArr(validIdx);
%     [~, maxIdx] = max(maxLevelActivations);
%     experts = projectNode(maxLevelNodes(maxIdx, 1:4), vocabularyDistributions, 'modal', options);
%     filterIds = round(((180/numel(options.filters)) * (0:(numel(options.filters)-1))) / (180/size(visFilters,3)))' + 1;
%     experts(:,1) = filterIds(experts(:,1));
%     modalImg = obtainPoE(experts, [], [], options.imageSize, visFilters, []);
%     
    %% Predict a category!
    idx = exportArr(:,4) == maxLevel;
    topLabels = exportArr(idx,1);
    uniqueTopLabels = unique(topLabels);
    combinedArr = cat(1, vocabulary{maxLevel}(uniqueTopLabels).categoryArr);
    combinedArr = mean(combinedArr, 1);
    [maxVal, predictedCategory] = max(combinedArr);
    % In case of more than 1 values, move on.
    if nnz(combinedArr == maxVal) > 1
         predictedCategory = datasample(find(combinedArr == maxVal), 1);
    end
    load([options.testInferenceFolder '/' categoryName '_' fileName '_test.mat']);
    
    %% Print realizations in the desired format for inte2D/3D integration.
    if exist([options.testInferenceFolder '/' categoryName '_' fileName '_test.mat'], 'file')
        save([options.testInferenceFolder '/' categoryName '_' fileName '_test.mat'], 'exportArr', 'activationArr', 'imgSize', 'predictedCategory', '-append');
    else
        save([options.testInferenceFolder '/' categoryName '_' fileName '_test.mat'], 'exportArr', 'activationArr', 'imgSize');
    end
end