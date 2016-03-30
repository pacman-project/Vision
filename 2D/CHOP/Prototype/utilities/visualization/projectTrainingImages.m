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
function [ ] = projectTrainingImages( fileList, vocabulary, mainGraph, levelItr, options)
    graphLevel = mainGraph{levelItr};
    imageIds = [graphLevel.imageId];
    realLabelIds = double(cat(1, graphLevel.realLabelId));
    precisePositions = double(cat(1, graphLevel.precisePosition));
    
    [rfSizes, visFilters, optimizedFilters, likelihoodLookupTable] = createOptimizationStructures(options, levelItr);
    
    %% Go through every image and find optimal version.
    for imgItr = 1:max(imageIds)
       % Obtain image-specific realizations.
       idx = imageIds == imgItr;
       
       % If no nodes are present, move on.
       if nnz(idx) == 0
            continue;
       end
       imagePrecisePositions = precisePositions(idx,:);

       % Create a folder if needed.
       [~, fileName, ~] = fileparts(fileList{imgItr});

       % Backproject from all possible levels.
       curExportArr = [realLabelIds(idx,:), imagePrecisePositions, repmat(levelItr, nnz(idx),1)];

       % Now, we get the top realizations and backproject to the original
       % image.
       optimizeImagination(curExportArr, vocabulary, options.imageSize, rfSizes, optimizedFilters, visFilters, 1, options.datasetName, likelihoodLookupTable, fileName);

%        % For visualization, overlay the original image with reconstructed nodes.
%        imwrite(muImg, [imgFolder '/' fileName '_level' num2str(levelItr) 'optimizedImagination.png']);
%        imwrite(orgImg, [imgFolder '/' fileName '_original.png']);
%        
    end
end
    