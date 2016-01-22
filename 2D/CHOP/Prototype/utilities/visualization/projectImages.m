%> Name: projectImages
%>
%> Description: Given the inference of parts in a dataset, this function
%> backprojects detected object models in real images. The backprojection
%> starts from the top level detections in every image, and continues all
%> the way to the bottom level. The same operation is performed both
%> training and test images.
%>
%> @param datasetName Name of the dataset to work on. 
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 19.11.2015
function [ ] = projectImages( datasetName)
    options = SetParameters(datasetName, false);
    outputFolder = [pwd '/output/' datasetName];
    numberOfFeatures = options.numberOfGaborFilters;
    radius = round(options.gaborFilterSize/2);
    
    filters = options.filterImages;
    
    %% Backproject training images.
    load([outputFolder '/export.mat'], 'exportArr', 'trainingFileNames', 'categoryNames', 'precisePositions');
    load([outputFolder '/vb.mat'], 'vocabulary');
    load([outputFolder '/mainGraph.mat'], 'mainGraph');
    realLabels = cellfun(@(x) [x.realLabelId], mainGraph, 'UniformOutput', false);
    realLabels = cat(2, realLabels{:})';
    exportArr(:,1) = realLabels;
    numberOfImages = max(exportArr(:,5));
    logLikelihoodArr = zeros(numberOfImages,1);
    for imgItr = 1:numberOfImages
       % Obtain image-specific realizations.
       idx = exportArr(:,5)==imgItr;
       imageExportArr = exportArr(idx,:);
       imagePrecisePositionArr = precisePositions(idx,:);
       
       % If there are no realizations, move on.
       if isempty(imageExportArr)
           continue;
       end
       maxLevel = max(imageExportArr(:,4));
       level1DataPoints = imageExportArr(imageExportArr(:,4) == 1, 1:3);
       
       [~, fileName, ~] = fileparts(trainingFileNames{imgItr});

       % Create a folder if needed.
       imgFolder = [pwd '/output/' datasetName '/reconstruction/train/' fileName];
       if exist(imgFolder, 'dir')
            rmdir(imgFolder, 's');
       end
       mkdir(imgFolder);

       % Get original image.
       img = imread(trainingFileNames{imgItr});
       
       % Backproject from all possible levels.
       for levelItr = 1:maxLevel
            idx = imageExportArr(:,4) == levelItr;
            curExportArr = imageExportArr(idx, :);
            curPrecisePositionArr = imagePrecisePositionArr(idx, :);
            curExportArr(:,2:3) = int32(round(curPrecisePositionArr));

            % Now, we get the top realizations and backproject to the original
            % image.
            level1Nodes = projectNode(curExportArr(:,1:4), vocabulary, 0, 'modal');
            level1NodesOrg = level1Nodes;
            if levelItr == maxLevel
                 % Calculate a mapping betweeen data points and reconstructed points.
                 [level1DataPoints, level1Nodes] = mapNodes(level1DataPoints, level1Nodes, radius);

                 % Finally, calculate log likelihood of the data. 
                 estimator = Estimator;
                 estimator = estimator.SetCoords(level1DataPoints(:,2:3), 0);
                 estimator = estimator.SetCoords(level1Nodes(:,2:3), 1);
                 estimator = estimator.SetLabels(level1DataPoints(:,1),0);
                 estimator = estimator.SetLabels(level1Nodes(:,1),1);
                 estimator = estimator.SetNumberOfFeatures(numberOfFeatures);
                 estimator = estimator.SetRadius(radius);
                 logLikelihoodArr(imgItr) = estimator.CalculateLogLikelihood();
            end
        
            % For visualization, overlay the original image with reconstructed nodes.
            projectionImg = overlayFeaturesWithImage(level1NodesOrg, img, filters, 1);
            imwrite(projectionImg, [pwd '/output/' datasetName '/reconstruction/train/' fileName '/' fileName '_1Projection' num2str(levelItr) '.png']);
            imwrite(img, [pwd '/output/' datasetName '/reconstruction/train/' fileName '/' fileName '_2Original.png']);
       end
    end
    trainLogLikelihood = sum(logLikelihoodArr);
    save([pwd '/output/' datasetName '/logLikelihood.mat'], 'trainLogLikelihood', 'logLikelihoodArr');
    
    %% Next, we do the same backprojection with test images.
    % First, we obtain the category information for every image.
    testFileNames = fuf([options.currentFolder '/input/' datasetName '/test/*.png'], 1, 'detail');
    categoryArrIdx = zeros(numel(testFileNames),1);
    for fileItr = 1:numel(testFileNames)
        categoryLabel = 1;
        fullName = testFileNames{fileItr};
        [~, fileName, ext] = fileparts(testFileNames{fileItr});
        strLength = numel([options.currentFolder '/input/' options.datasetName '/test/']);
        fileNameLength = numel(ext) + numel(fileName) + 1; % 1 for '/' (folder delimiter)
        if numel(fullName) >= strLength + fileNameLength
            categoryStr = fullName(1, (strLength+1):(end-fileNameLength));
            categoryStrSepIdx = strfind(categoryStr, '/');
            if ~isempty(categoryStrSepIdx)
                categoryStr = categoryStr(:, 1:(categoryStrSepIdx(1)-1));
            end
            chosenCategoryArr = cellfun(@(x) strcmp(x, categoryStr), categoryNames);
            categoryLabel = find(chosenCategoryArr, 1, 'first');
        end
        if ~isempty(categoryLabel)
            categoryArrIdx(fileItr) = categoryLabel;
        end
    end

    % Next, we find the output files and do backprojections.
    numberOfImages = numel(testFileNames);
    testLogLikelihoodArr = zeros(numberOfImages,1);
    if ~exist([pwd '/output/' datasetName '/reconstruction/test'], 'dir')
        mkdir([pwd '/output/' datasetName '/reconstruction/test']);
    end
    for imgItr = 1:numberOfImages
       [~, testFileName, ~] = fileparts(testFileNames{imgItr});
       outputFile = [options.testInferenceFolder '/' categoryNames{categoryArrIdx(imgItr)} '_' testFileName '_test.mat'];
       load(outputFile);
       exportArr(:,1) = realLabelIds;
       % Obtain image-specific realizations.
        % Assign real labels.
       imageExportArr = exportArr;
       
       % If there are no realizations, move on.
       if isempty(imageExportArr)
           continue;
       end
       
      % Create a folder if needed.
       imgFolder = [pwd '/output/' datasetName '/reconstruction/test/' testFileName];
       if exist(imgFolder, 'dir')
            rmdir(imgFolder, 's');
       end
       mkdir(imgFolder);
       
       % Backproject from all levels.
       maxLevel = max(imageExportArr(:,4));
       level1DataPoints = imageExportArr(imageExportArr(:,4) == 1, 1:3);
       

       for levelItr = 1:maxLevel
            idx = imageExportArr(:,4) == levelItr;
            curExportArr = imageExportArr(idx, :);
            curExportArr(:,2:3) = int32(round(precisePositions(idx, :)));

            % Get original image.
            [~, fileName, ~] = fileparts(testFileNames{imgItr});
            img = imread(testFileNames{imgItr});

            % Now, we get the top realizations and backproject to the original
            % image.
            level1Nodes = projectNode(curExportArr(:,1:4), vocabulary, 0, 'modal');
            level1NodesOrg = level1Nodes;

            if levelItr == maxLevel
                 % Calculate a mapping betweeen data points and reconstructed points.
                 [level1DataPoints, level1Nodes] = mapNodes(level1DataPoints, level1Nodes, radius);
                 
                 % Finally, calculate log likelihood of the data. 
                 estimator = Estimator;
                 estimator = estimator.SetCoords(level1DataPoints(:,2:3), 0);
                 estimator = estimator.SetCoords(level1Nodes(:,2:3), 1);
                 estimator = estimator.SetLabels(level1DataPoints(:,1),0);
                 estimator = estimator.SetLabels(level1Nodes(:,1),1);
                 estimator = estimator.SetNumberOfFeatures(numberOfFeatures);
                 estimator = estimator.SetRadius(radius);
                 testLogLikelihoodArr(imgItr) = estimator.CalculateLogLikelihood();
            end
            
            % For visualization, overlay the original image with reconstructed nodes.
            projectionImg = overlayFeaturesWithImage(level1NodesOrg, img, filters, 1);
            imwrite(projectionImg, [imgFolder '/' fileName '_1Projection' num2str(levelItr) '.png']);
            imwrite(img, [imgFolder '/' fileName '_2Original.png']);
       end
    end
    testLogLikelihood = sum(testLogLikelihoodArr);
    save([pwd '/output/' datasetName '/logLikelihood.mat'], 'testLogLikelihood', 'testLogLikelihoodArr', '-append');
end

function [level1DataPoints, level1NodesNew] = mapNodes(level1DataPoints, level1Nodes, radius)

    numberOfDataPoints = size(level1DataPoints,1);
    mappingIdx = zeros(numberOfDataPoints,1);
    distances = zeros(numberOfDataPoints,1);
    
    D = pdist2(double(level1DataPoints(:,2:3)), double(level1Nodes(:,2:3)));
    for dataPointItr = 1:numberOfDataPoints
        [Y, dataId] = min(D);
        [distance, reconId] = min(Y);
        dataId = dataId(reconId);
        if distance > radius
            mappingIdx(dataId) = 0;
        else
            mappingIdx(dataId) = reconId;
            distances(dataId) = distance;
        end
        % dataId is the id of the reconstructed node, while reconId is the data
        % point assigned. 
        D(dataId,:) = inf;
        D(:, reconId) = inf;
    end
    assignedNodeIdx = mappingIdx>0;
    assignedReconNodeIdx = mappingIdx(mappingIdx>0);
    level1NodesNew = zeros(size(level1DataPoints), 'int32');
    level1NodesNew(assignedNodeIdx,:) = level1Nodes(assignedReconNodeIdx,:);
    level1DataPoints = uint16(level1DataPoints);
    level1NodesNew = uint16(level1NodesNew);
end