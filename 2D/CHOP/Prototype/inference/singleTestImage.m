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
function [] = singleTestImage(testFileName, vocabulary, vocabularyDistributions, allModes, modeProbs, categoryName, options)
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
    [ nodes, ~, nodeActivations, ~, ~, ~ ] = getNodes( img, [], options );
    if isempty(nodes)
        return;
    end
    % Save smoothed image.
    % Assign nodes their image ids.
    nodes = int32(cell2mat(nodes));
    [exportArr, activationArr] = inferSubs(img, vocabulary, vocabularyDistributions, allModes, modeProbs, nodes, nodeActivations, options); %#ok<ASGLU,NASGU>
    maxLevel = max(exportArr(:,4));
    
    %% Project stuff from top layer.    % Create data structures required for optimization.
    if ~exist([pwd '/filters/optimizationFilters.mat'], 'file')
         [rfSizes, visFilters, optimizedFilters, likelihoodLookupTable] = createOptimizationStructures(options, levelItr, true);
         save([pwd '/filters/optimizationFilters.mat'], 'rfSizes', 'visFilters', 'optimizedFilters', 'likelihoodLookupTable');
    else
         load([pwd '/filters/optimizationFilters.mat'], 'visFilters');
    end
    validIdx = exportArr(:,4) == maxLevel;
    maxLevelNodes = exportArr(validIdx,:);
    maxLevelActivations = activationArr(validIdx);
    [~, maxIdx] = max(maxLevelActivations);
    experts = projectNode(maxLevelNodes(maxIdx, 1:4), vocabularyDistributions, 'modal');
    filterIds = round(((180/numel(options.filters)) * (0:(numel(options.filters)-1))) / (180/size(visFilters,3)))' + 1;
    experts(:,1) = filterIds(experts(:,1));
    modalImg = obtainPoE(experts, [], [], options.imageSize, visFilters, []);
    
    %% Print realizations in the desired format for inte2D/3D integration.
    if exist([options.testInferenceFolder '/' categoryName '_' fileName '_test.mat'], 'file')
        save([options.testInferenceFolder '/' categoryName '_' fileName '_test.mat'], 'exportArr', 'activationArr', 'imgSize', 'precisePositions', 'realLabelIds', '-append');
    else
        save([options.testInferenceFolder '/' categoryName '_' fileName '_test.mat'], 'exportArr', 'activationArr', 'imgSize', 'precisePositions', 'realLabelIds');
    end
end