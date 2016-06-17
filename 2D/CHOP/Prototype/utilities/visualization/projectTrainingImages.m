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
    % Create data structures required for optimization.
    if ~exist([pwd '/filters/optimizationFilters.mat'], 'file')
         [rfSizes, visFilters, optimizedFilters, likelihoodLookupTable] = createOptimizationStructures(options, levelItr, true);
         save([pwd '/filters/optimizationFilters.mat'], 'rfSizes', 'visFilters', 'optimizedFilters', 'likelihoodLookupTable');
    else
         load([pwd '/filters/optimizationFilters.mat'], 'rfSizes', 'visFilters', 'optimizedFilters', 'likelihoodLookupTable');
    end
    levelExportArr = exportArr(exportArr(:,4) == levelItr, :);
    
    %% Create optimization options.
     optimizationOptions.stopVal = 0.1;
     optimizationOptions.maxSteps = 1;
     optimizationOptions.minOptimizationLayer = 3;
     optimizationOptions.minLikelihoodChange = 0.000001;
     optimizationOptions.movesPerChild = 3;
     optimizationOptions.maxMoves = 10;
     optimizationOptions.minMoves = 5;
     optimizationOptions.likelihoodChangeThr =  1.000001;
     optimizationOptions.moveFlags = [1,1,0]; % 1 for position moves, 2 is for or moves, 3 for rotation moves.
     % If poeFlag is true, we are searching for pixel level agreement.
     optimizationOptions.poeFlag = true;
     % Position flags are the strings that keep things in place.
     optimizationOptions.positionFlag = true;
    
    %% Go through every image and find optimal version.
    parfor imgItr = 1:max(levelExportArr(:,5))
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
       [ experts, ~, ~, ~, ~ ] = projectNode(curExportArr, vocabularyDistributions, 'modal', options);
       
       % If rotation flag is on, we'll consider a wider range of filters.
       numberOfFilters = size(visFilters,3);
       numberOfRealFilters = numel(vocabularyDistributions{1});
       if numberOfFilters > numel(vocabularyDistributions{1})
             filterIds = round(((180/numberOfRealFilters) * (0:(numberOfRealFilters-1))) / (180/numberOfFilters))' + 1;
             experts(:,1) = filterIds(experts(:,1));
       end
       
       % Visualize the experts.
       [refModalImg, ~, ~] = obtainPoE(experts, [], [], options.imageSize, visFilters, []);
       
%       refModalImg = optimizeImagination(curExportArr, vocabularyDistributions, options.imageSize, rfSizes, optimizedFilters, visFilters, 1, options.datasetName, likelihoodLookupTable, options, optimizationOptions, fileName);

%        % For visualization, overlay the original image with reconstructed nodes.
        imwrite(refModalImg, [options.outputFolder '/reconstruction/train/' fileName '/' fileName '_level' num2str(levelItr) 'imagination.png']);
%        imwrite(orgImg, [imgFolder '/' fileName '_original.png']);
%        
    end
end
    