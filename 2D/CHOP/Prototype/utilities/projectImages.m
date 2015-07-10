function [ ] = projectImages( datasetName )
    options = SetParameters(datasetName, false);
    outputFolder = [pwd '/output/' datasetName];
    if strcmp(options.filterType, 'gabor')
       inhibitionRadius = options.gabor.inhibitionRadius; 
    else
       inhibitionRadius = options.auto.inhibitionRadius;
    end
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
    for imgItr = 1:numberOfImages
       % Obtain image-specific realizations.
       imageExportArr = exportArr(exportArr(:,5)==imgItr,:);
       
       % If there are no realizations, move on.
       if isempty(imageExportArr)
           continue;
       end
       
       maxLevel = max(imageExportArr(:,4));
       imageExportArr = imageExportArr(imageExportArr(:,4) == maxLevel, :);
       
       % Get original image.
       [~, fileName, ~] = fileparts(trainingFileNames{imgItr});
       img = imread(trainingFileNames{imgItr});
       
       % Now, we get the top realizations and backproject to the original
       % image.
       level1Nodes = projectNode(imageExportArr(:,1:4), vocabulary, inhibitionRadius);
        
       % Finally, overlay the original image with reconstructed nodes.
       img = overlayFeaturesWithImage(level1Nodes, img, filters);
       imwrite(img, [pwd '/output/' datasetName '/reconstruction/train/' fileName '_backProjection.png']);
    end
end

