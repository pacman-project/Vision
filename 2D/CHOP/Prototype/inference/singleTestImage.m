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
function [] = singleTestImage(testFileName, vocabulary, distanceMatrices, optimalThresholds, options)
    %% Get the first level nodes.
    % First, downsample the image if it is too big.
    img = imread(testFileName);
    [~, fileName, ~] = fileparts(testFileName);
    if options.debug
        display(['Processing ' fileName '.']);
    end
    % resize image if necessary.
    if max(size(img)) > options.maxImageDim
       img = imresize(img, options.maxImageDim/max(size(img)), 'bilinear'); 
    end
%    imwrite(img, [options.processedFolder '/' fileName '.png']);

    %% Form the first level nodes.
    [cellNodes, ~] = getNodes(img, [], options);
%    imwrite(smoothedImg, [options.smoothedFolder '/' fileName '.png']);
    if isempty(cellNodes)
        return;
    end
    % Save smoothed image.
    % Assign nodes their image ids.
    nodes = zeros(size(cellNodes,1), 3, 'int32');
    nodes(:,1:3) = cell2mat(cellNodes);
    [exportArr, confidenceArr] = inferSubs(vocabulary, nodes, distanceMatrices, optimalThresholds, options);
    
    %% Print realizations in the desired format for inte2D/3D integration.
    if exist([options.testInferenceFolder '/' fileName '_test.mat'], 'file')
        save([options.testInferenceFolder '/' fileName '_test.mat'], 'exportArr', 'confidenceArr', '-append');
    else
        save([options.testInferenceFolder '/' fileName '_test.mat'], 'exportArr', 'confidenceArr');
    end
end