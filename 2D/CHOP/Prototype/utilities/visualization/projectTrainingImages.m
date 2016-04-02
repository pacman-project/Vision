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
function [ ] = projectTrainingImages( fileList, vocabularyDistributions, exportArr, levelItr, options)

    [rfSizes, visFilters, optimizedFilters, likelihoodLookupTable] = createOptimizationStructures(options, levelItr, true);
    levelExportArr = exportArr(exportArr(:,4) == levelItr, :);
    
    %% Go through every image and find optimal version.
    for imgItr = 10:max(levelExportArr(:,5))
       % Create a folder if needed.
       [~, fileName, ~] = fileparts(fileList{imgItr});

       % Backproject from all possible levels.
       curExportArr = levelExportArr(levelExportArr(:,5) == imgItr,:);
       
       % If no nodes are present, move on.
       if isempty(curExportArr)
            continue;
       end

       % Now, we get the top realizations and backproject to the original
       % image.
       optimizeImagination(curExportArr, vocabularyDistributions, options.imageSize, rfSizes, optimizedFilters, visFilters, 1, options.datasetName, likelihoodLookupTable, fileName);

%        % For visualization, overlay the original image with reconstructed nodes.
%        imwrite(muImg, [imgFolder '/' fileName '_level' num2str(levelItr) 'optimizedImagination.png']);
%        imwrite(orgImg, [imgFolder '/' fileName '_original.png']);
%        
    end
end
    