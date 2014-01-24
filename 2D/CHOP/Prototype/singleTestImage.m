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
function [] = singleTestImage(testFileName, options, currentPath)
    % Here, we will run the inference process by compressing the test
    % images' graphs with the compositions in the vocabulary.
    % Allocate space for current graph level.
    load([currentPath '/output/' options.datasetName '/' options.datasetName '_vb.mat'], 'vocabulary', 'modes');
    mainGraph = cell(options.maxLevels,1);
    
    %% Get the first level nodes.
    % First, downsample the image if it is too big.
    img = imread(testFileName);
    [~, fileName, ~] = fileparts(testFileName);
    if options.debug
        display(['Processing ' fileName '.']);
    end
    if max(size(img)) > options.maxImageDim
       img = imresize(img, options.maxImageDim/max(size(img)), 'bilinear'); 
    end
    imwrite(img, [options.processedFolder '/' fileName '.png']);

    %% Form the first level nodes.
    nodes = getNodes(img, options);
    % Assign nodes their image ids.
    imageIds = ones(size(nodes,1), 1);
    nodes(:, 3) = mat2cell(imageIds, ones(size(imageIds)));
    graphLevel(size(nodes,1)) = struct('labelId', [], 'imageId', [], 'position', [], 'children', [], 'parents', [], 'adjInfo', [], 'leafNodes', []);

    %% Get edges depending on the property to be embedded in the graph.
    [~, edges, leafNodeAdjArr] = extractEdges(nodes, [], [], options, 1, options.datasetName, modes);
    %% Fill the basic info in this scene graph level.
    for instanceItr = 1:size(nodes,1)
       graphLevel(instanceItr).labelId = nodes{instanceItr,1};
       graphLevel(instanceItr).position = fix(nodes{instanceItr,2});
       graphLevel(instanceItr).imageId = imageIds(instanceItr);
       graphLevel(instanceItr).leafNodes = instanceItr;
       nodeEdges = ismember(edges(:,1:2), instanceItr);

       % get non-zero rows of edges
       selfEdges = edges(nodeEdges(:,1) | nodeEdges(:,2),:);
       if numel(selfEdges)>0
            graphLevel(instanceItr).adjInfo = selfEdges;
       end
    end
    
    %% Assign first level of the scene graph.
    mainGraph(1) = {graphLevel};
    firstLevel = graphLevel;
    
    % Create folder to put the output structures in.
    if exist([options.testGraphFolder '/' fileName], 'dir')
       rmdir([options.testGraphFolder '/' fileName], 's');
    end
    mkdir([options.testGraphFolder '/' fileName]);
    
    %% Iteratively process each level to parse the object.
    for levelItr = 2:numel(vocabulary)
        graphFileName = [options.testGraphFolder '/' fileName '/level' num2str(levelItr-1) '.g'];
        fp = fopen(graphFileName, 'w');
        
        %% Visualize the test images with previous layer's subs.
        if options.debug
            visualizeImages({testFileName}, mainGraph, levelItr-1, options, options.datasetName, 'test' );
        end
        
        %% Print the graph to the input file.
        printGraphToFile(fp, nodes(:,1), edges, true);
        fclose(fp);
        
        %% Here, we run SUBDUE over the input graph(s) to find pre-defined compositions within the graph.
        % Each pre-defined sub is searched separately.
        if options.debug
           display(['Working on level ' num2str(levelItr) '.']);
        end
        newLevel = collectInstances(vocabulary{levelItr}, graphFileName, options);
        
        %% Assign positions, image ids, and leaf nodes. 
        % If no new subs have been found, finish processing.
        if isempty(newLevel)
           mainGraph = mainGraph(1:(levelItr-1),:);
           break; 
        end
        
        % Assign positions, image ids and leaf nodes.
        previousLevel = mainGraph{levelItr-1};
        for newNodeItr = 1:numel(newLevel)
            leafNodes = [];
            position = [0,0];
            newLevel(newNodeItr).imageId = previousLevel(newLevel(newNodeItr).children(1)).imageId;
            for childItr = 1:numel(newLevel(newNodeItr).children)
               leafNodes = [leafNodes, previousLevel(newLevel(newNodeItr).children(childItr)).leafNodes]; 
            end
            leafNodes = unique(leafNodes);
            for leafNodeItr = 1:numel(leafNodes)
               position = position + firstLevel(leafNodes(leafNodeItr)).position;
            end
            newLevel(newNodeItr).leafNodes = leafNodes;
            newLevel(newNodeItr).position = round(position/numel(leafNodes));
        end
        
        %% Apply local inhibition.
        [newLevel] = applyLocalInhibition(newLevel, options, levelItr);
        
        %% If new level is empty, break.
        if isempty(newLevel)
            break;
        end
        
        %% Create new graph to be fed to for knowledge discovery in the next level.
        numberOfNodes = numel(newLevel);
        nodes = cell(numel(newLevel), 3);
        for nodeItr = 1:numberOfNodes
           nodes(nodeItr,:) = {newLevel(nodeItr).labelId, newLevel(nodeItr).position, newLevel(nodeItr).imageId};
        end
        
        %% Create parent relationships.
        mainGraph = mergeIntoGraph(mainGraph, newLevel, levelItr, 1);
        
        %% Extract the edges to form the new graph.
        if numel(modes)<levelItr
           edges=[];
        else
            [~, edges] = extractEdges(nodes, mainGraph, leafNodeAdjArr, options, levelItr, options.datasetName, modes);
        end
        
        %% If the edges are not empty, fill the edge information in current level.
        if ~isempty(edges)
            newLevel = mainGraph{levelItr};
            for instanceItr = 1:size(nodes,1)
               nodeEdges = ismember(edges(:,1:2), instanceItr);

               % get non-zero rows of edges
               selfEdges = edges(nodeEdges(:,1) | nodeEdges(:,2),:);
               if numel(selfEdges)>0
                    newLevel(instanceItr).adjInfo = selfEdges;
               end
            end
            
            mainGraph{levelItr} = newLevel;
        end
        
        %% Visualize the test images with previous layer's subs.
        if options.debug && levelItr == numel(vocabulary)
            visualizeImages( {testFileName}, mainGraph, levelItr, options, options.datasetName, 'test' );
        end
    end
    filledIdx = find(~(cellfun('isempty', mainGraph)), 1, 'first');
    mainGraph = mainGraph(1:filledIdx);
    save([currentPath '/output/' options.datasetName '/' options.datasetName '_test.mat'], 'mainGraph');
end

