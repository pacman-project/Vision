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
function [] = singleTestImage(testFileName, vocabulary, vocabularyDistributions, categoryName, options)
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
    [exportArr, activationArr] = inferSubs(vocabulary, vocabularyDistributions, nodes, nodeActivations, options); %#ok<ASGLU,NASGU>
    
    %% Print realizations in the desired format for inte2D/3D integration.
    if exist([options.testInferenceFolder '/' categoryName '_' fileName '_test.mat'], 'file')
        save([options.testInferenceFolder '/' categoryName '_' fileName '_test.mat'], 'exportArr', 'activationArr', 'imgSize', 'precisePositions', 'realLabelIds', '-append');
    else
        save([options.testInferenceFolder '/' categoryName '_' fileName '_test.mat'], 'exportArr', 'activationArr', 'imgSize', 'precisePositions', 'realLabelIds');
    end
end