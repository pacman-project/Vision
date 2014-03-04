%> Name: singleTestImages
%>
%> Description: Process each image, run discovery with vocabulary at each
%> level, and find the class each image belongs to. The difference from
%testImages is that each image is processed separately.
%>
%> @param testFileImages The test image names to work on.
%> @param options Program options.
%>
%> @retval classes The classification results, the class of each given
%> image.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 19.12.2013
%> 
function [] = singleTestImage(testFileName, options)
    % Here, we will run the inference process by compressing the test
    % images' graphs with the compositions in the vocabulary.
    % Allocate space for current graph level.
    load([options.currentFolder '/output/' options.datasetName '/' options.datasetName '_vb.mat'], 'vocabulary', 'modes', 'highLevelModes');
    
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
    nodes = getNodes(img, [], options);
    % Assign nodes their image ids.
    imageIds = ones(size(nodes,1), 1);
    nodes(:, 3) = mat2cell(imageIds, ones(size(imageIds)));
    leafNodes = nodes;
    
    % Generate first level object graph.
    [~, graphLevel] = generateLevels(nodes, options);

    %% Step 2.1: Get first-level object graph edges.
    mainGraph = {graphLevel};
    
    %% Get edges depending on the property to be embedded in the graph.
    [~, ~, mainGraph] = extractEdges(mainGraph, options, 1, modes, highLevelModes);
    
    %% Visualize level 1 test image.
    if options.debug
        visualizeImages( {testFileName}, mainGraph{1}, leafNodes, 1, options, 'test' );
    end
    
    %% Iteratively process each level to parse the object.
    for levelItr = 2:numel(vocabulary)
        %% Here, we run SUBDUE over the input graph(s) to find pre-defined compositions within the graph.
        % Each pre-defined sub is searched separately.
        if options.debug
           display(['Working on level ' num2str(levelItr) '.']);
        end
        newLevel = collectInstances(vocabulary{levelItr}, mainGraph{levelItr-1}, [], options, levelItr);
        
        %% Assign positions, image ids, and leaf nodes. 
        % If no new subs have been found, finish processing.
        if isempty(newLevel)
           break; 
        end
        
        %% Fill in children, position info and sort nodes based on image id.
        previousLevel = mainGraph{levelItr-1};
        newLevel = fillBasicInfo(previousLevel, newLevel, leafNodes);
        
        %% Apply local inhibition.
        [newLevel] = applyTestInhibition(newLevel, options, levelItr);
%        [newLevel] = applyLocalInhibition(newLevel, options, levelItr);(vocabLevel, graphLevel, currentModes, options, levelItr)
        
        %% If new level is empty, break.
        if isempty(newLevel)
            if options.debug
                visualizeImages( {testFileName}, mainGraph{levelItr}, leafNodes, levelItr, options, 'test' );
            end
            break;
        end
        
        %% Create parent relationships.
        mainGraph = mergeIntoGraph(mainGraph, newLevel, [], levelItr, 1);
        
        %% Extract the edges to form the new graph.
        [~, ~, mainGraph] = extractEdges(mainGraph, options, levelItr, modes, highLevelModes);
        
        %% Visualize the test images with previous layer's subs.
        if options.debug
            visualizeImages( {testFileName}, mainGraph{levelItr}, leafNodes, levelItr, options, 'test' );
        end
    end
    save([options.testInferenceFolder '/' fileName '_test.mat'], 'mainGraph');
end

