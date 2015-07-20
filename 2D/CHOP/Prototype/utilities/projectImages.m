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
    load([outputFolder '/export.mat'], 'exportArr', 'trainingFileNames');
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
       img = overlayFeaturesWithImage(level1Nodes, img, filters);
       imwrite(img, [pwd '/output/' datasetName '/reconstruction/train/' fileName '_backProjection.png']);
    end
    logLikelihood = sum(logLikelihoodArr);
    save([pwd '/output/' datasetName '/logLikelihood.mat'], 'logLikelihood', 'logLikelihoodArr');
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