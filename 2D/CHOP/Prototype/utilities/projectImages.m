function [ ] = projectImages( datasetName)
    options = SetParameters(datasetName, false);
    outputFolder = [pwd '/output/' datasetName];
    if strcmp(options.filterType, 'gabor')
       inhibitionRadius = options.gabor.inhibitionRadius; 
    else
       inhibitionRadius = options.auto.inhibitionRadius;
    end
    numberOfFeatures = options.numberOfGaborFilters;
    radius = round(options.gaborFilterSize/2);
    
    filters = options.filters;
    
    %% Backproject training images.
    load([outputFolder '/export.mat'], 'exportArr', 'trainingFileNames', 'categoryNames');
    load([outputFolder '/vb.mat'], 'vocabulary');
    filters = cellfun(@(x) uint8(round(255 * (x - min(min(min(x)))) / (max(max(max(x))) - min(min(min(x)))))), filters, 'UniformOutput', false);
    img = imread(trainingFileNames{1});
    bandCount = size(img,3);
    if size(filters{1}, 3) < bandCount
        filters = cellfun(@(x) repmat(x, [1 1 bandCount]), filters, 'UniformOutput', false);
    end
    
    numberOfImages = max(exportArr(:,5));
    logLikelihoodArr = zeros(numberOfImages,1);
    for imgItr = 1:numberOfImages
       % Obtain image-specific realizations.
       imageExportArr = exportArr(exportArr(:,5)==imgItr,:);
       
       % If there are no realizations, move on.
       if isempty(imageExportArr)
           continue;
       end
       
       maxLevel = max(imageExportArr(:,4));
       level1DataPoints = imageExportArr(imageExportArr(:,4) == 1, 1:3);
       imageExportArr = imageExportArr(imageExportArr(:,4) == maxLevel, :);
       
       % Get original image.
       [~, fileName, ~] = fileparts(trainingFileNames{imgItr});
       img = imread(trainingFileNames{imgItr});
       
       % Now, we get the top realizations and backproject to the original
       % image.
       level1Nodes = projectNode(imageExportArr(:,1:4), vocabulary, inhibitionRadius);
       level1NodesOrg = level1Nodes;
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
        
       % For visualization, overlay the original image with reconstructed nodes.
       img = overlayFeaturesWithImage(level1NodesOrg, img, filters);
       imwrite(img, [pwd '/output/' datasetName '/reconstruction/train/' fileName '_backProjection.png']);
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
       load(outputFile, 'exportArr');
       % Obtain image-specific realizations.
       imageExportArr = exportArr;
       
       % If there are no realizations, move on.
       if isempty(imageExportArr)
           continue;
       end
       
       maxLevel = max(imageExportArr(:,4));
       level1DataPoints = imageExportArr(imageExportArr(:,4) == 1, 1:3);
       imageExportArr = imageExportArr(imageExportArr(:,4) == maxLevel, :);
       
       % Get original image.
       [~, fileName, ~] = fileparts(testFileNames{imgItr});
       img = imread(testFileNames{imgItr});
       
       % Now, we get the top realizations and backproject to the original
       % image.
       level1Nodes = projectNode(imageExportArr(:,1:4), vocabulary, inhibitionRadius);
       level1NodesOrg = level1Nodes;
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
        
       % For visualization, overlay the original image with reconstructed nodes.
       img = overlayFeaturesWithImage(level1NodesOrg, img, filters);
       imwrite(img, [pwd '/output/' datasetName '/reconstruction/test/' fileName '_backProjection.png']);
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