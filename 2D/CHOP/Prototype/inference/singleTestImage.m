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
function [] = singleTestImage(testFileName, vocabulary, distanceMatrices, graphLevelIndices, options)
    %% Get the first level nodes.
    % First, downsample the image if it is too big.
    img = imread(testFileName);
    threshold = options.subdue.threshold;
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
    [cellNodes, smoothedImg] = getNodes(img, [], options);
    imwrite(smoothedImg, [options.smoothedFolder '/' fileName '.png']);
    if isempty(cellNodes)
        return;
    end
    % Save smoothed image.
    % Assign nodes their image ids.
    imageIds = ones(size(cellNodes,1), 1, 'int32');
    nodes = zeros(size(cellNodes,1), 4, 'int32');
    nodes(:,1:3) = cell2mat(cellNodes);
    nodes(:, 4) = imageIds;
    leafNodes = nodes;
    
    %% Generate first level object graph.
    [~, newLevel] = generateLevels(nodes, ones(size(nodes,1),1), options);

    %% Step 2.1: Get first-level object graph edges.
    mainGraph = {newLevel};
    
    %% Get edges depending on the property to be embedded in the graph.
    [mainGraph] = extractEdges(mainGraph, options, 1);
    
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
        newLevel = collectInstances(vocabulary{levelItr}, mainGraph{levelItr-1}, distanceMatrices{levelItr-1}, options, threshold, levelItr);
        
        %% Assign positions, image ids, and leaf nodes. 
        % If no new subs have been found, finish processing.
        if isempty(newLevel)
           break; 
        end
        
        %% Fill in children, position info and sort nodes based on image id.
        previousLevel = mainGraph{levelItr-1};
        newLevel = fillBasicInfo(previousLevel, newLevel, leafNodes, options.numberOfThreads);
        
        %% Apply local inhibition.
        if options.fastInference
            display('........ Applying inhibition.');
            [newLevel] = applyTestInhibition(newLevel, options, levelItr);
            display(['........ Inhibition applied with novelty thr: ' num2str(options.noveltyThr) ' and edge novelty thr: ' num2str(options.edgeNoveltyThr) '.']);
            display(['........ Remaining: ' num2str(numel(newLevel)) ' realizations belonging to ' num2str(numel(unique([newLevel.labelId]))) ' compositions.']);
        end
        
        %% Assign new labels to newLevel, and continue.    
        labelIds = [newLevel.labelId];
        graphLabelAssgnArr = graphLevelIndices{levelItr};
        newLabelIds = num2cell(graphLabelAssgnArr(labelIds)');
        [newLevel.labelId] = deal(newLabelIds{:});

        % Rearrange graph level so it is sorted by image id.
        arrayToSort = [[newLevel.imageId]', [newLevel.labelId]'];
        [~, sortedIdx] = sortrows(arrayToSort);
        newLevel = newLevel(sortedIdx);
        
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
        if levelItr < numberOfLevels
            [mainGraph] = extractEdges(mainGraph, options, levelItr);
        end
        
        %% Visualize the test images with previous layer's subs.
        if options.debug
            visualizeImages( {testFileName}, vocabulary{levelItr}, mainGraph{levelItr}, leafNodes, levelItr, options, 'test' );
        end
    end
    
    %% Process mainGraph to export realizations in the desired format for inte2D/3D integration.
    exportArr = exportRealizations(mainGraph); %#ok<NASGU>
    if exist([options.testInferenceFolder '/' fileName '_test.mat'], 'file')
        save([options.testInferenceFolder '/' fileName '_test.mat'], 'exportArr', '-append');
    else
        save([options.testInferenceFolder '/' fileName '_test.mat'], 'exportArr');
    end
end