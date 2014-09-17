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
function [totalInferenceTime] = singleTestImage(testFileName, vocabulary, redundantVocabulary, modes, options)
    totalInferenceTime = 0;
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
    imwrite(img, [options.processedFolder '/' fileName '.png']);

    %% Form the first level nodes.
    [nodes, smoothedImg] = getNodes(img, [], options);
    % Save smoothed image.
    imwrite(smoothedImg, [options.smoothedFolder '/' fileName '.png']);
    % Assign nodes their image ids.
    imageIds = ones(size(nodes,1), 1);
    nodes(:, 3) = mat2cell(imageIds, ones(size(imageIds)));
    leafNodes = nodes;
    
    if isempty(nodes)
        return;
    end
    
    %% Generate first level object graph.
    [~, graphLevel] = generateLevels(nodes, ones(size(nodes,1),1), options);

    %% Step 2.1: Get first-level object graph edges.
    mainGraph = {graphLevel};
    
    %% Get edges depending on the property to be embedded in the graph.
    [~, mainGraph] = extractEdges(mainGraph, options, 1, modes);
    
    %% Visualize level 1 test image.
    if options.debug
        visualizeImages( {testFileName}, vocabulary{1}, mainGraph{1}, leafNodes, 1, options, 'test' );
    end
    
    %% Iteratively process each level to parse the object.
    numberOfLevels = min(options.maxInferenceLevels, numel(vocabulary));
    for levelItr = 2:numberOfLevels
        %% Here, we run SUBDUE over the input graph(s) to find pre-defined compositions within the graph.
        % Each pre-defined sub is searched separately.
        if options.debug
           display(['Working on level ' num2str(levelItr) '.']);
        end
        startTime = tic;
        newLevel = collectInstances(vocabulary{levelItr}, redundantVocabulary{levelItr}, mainGraph{levelItr-1}, options, levelItr);
        duration = toc(startTime);
        totalInferenceTime = totalInferenceTime + duration;
        
        %% Assign positions, image ids, and leaf nodes. 
        % If no new subs have been found, finish processing.
        if isempty(newLevel)
           break; 
        end
        
        %% Fill in children, position info and sort nodes based on image id.
        previousLevel = mainGraph{levelItr-1};
        newLevel = fillBasicInfo(previousLevel, newLevel, leafNodes);
        
        %% Apply local inhibition.
        if options.fastInference
            display('........ Applying inhibition.');
            [newLevel] = applyTestInhibition(newLevel, options, levelItr);
            display(['........ Inhibition applied with novelty thr: ' num2str(options.noveltyThr) ' and edge novelty thr: ' num2str(options.edgeNoveltyThr) '.']);
            display(['........ Remaining: ' num2str(numel(newLevel)) ' realizations belonging to ' num2str(numel(unique([newLevel.labelId]))) ' compositions.']);
        end
        
        %% If new level is empty, break.
        if isempty(newLevel)
            if options.debug
                visualizeImages( {testFileName}, vocabulary{levelItr}, mainGraph{levelItr}, leafNodes, levelItr, options, 'test' );
            end
            break;
        end
        
        %% Create parent relationships.
        mainGraph = mergeIntoGraph(mainGraph, newLevel, [], levelItr, 1);
        
        %% Extract the edges to form the new graph.
        if levelItr ~= numberOfLevels
            [~, mainGraph] = extractEdges(mainGraph, options, levelItr, modes);
        end
        
        %% Visualize the test images with previous layer's subs.
        if options.debug
            visualizeImages( {testFileName}, vocabulary{levelItr}, mainGraph{levelItr}, leafNodes, levelItr, options, 'test' );
        end
    end
    
    %% Process mainGraph to export realizations in the desired format for inte2D/3D integration.
    exportArr = exportRealizations(mainGraph); %#ok<NASGU>
    if exist([options.testInferenceFolder '/' fileName '_test.mat'], 'file')
        save([options.testInferenceFolder '/' fileName '_test.mat'], 'exportArr');
    else
        save([options.testInferenceFolder '/' fileName '_test.mat'], 'exportArr');
    end
end