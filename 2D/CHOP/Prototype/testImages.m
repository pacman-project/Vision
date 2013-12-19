%> Name: testImages
%>
%> Description: Process each image, run discovery with vocabulary at each
%> level, and find the class each image belongs to.
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
%> Ver 1.0 on 15.12.2013
function [ classes ] = testImages(testFileNames, modes, leafNodeAdjArr, options, currentPath)
    numberOfTestImages = size(testFileNames,1);
    classes = zeros(numberOfTestImages,1);
    
    % Here, we will run the inference process by compressing the test
    % images' graphs with the compositions in the vocabulary.
    %% Step 1: Extract nodes of each image for the first level of the hierarchy.
    nodeCounter = 0;
    allNodes = cell(options.maxNumberOfFeatures,3);
    % Allocate space for current graph level.
    load([currentPath '/output/' options.datasetName '_vb.mat'], 'vocabulary', 'mainGraph', 'modes');
    mainGraph = cell(options.maxLevels,1);
    
    %% Get the first level nodes.
    for fileItr = 1:size(testFileNames,1) 
        %% First, downsample the image if it is too big.
        img = imread(testFileNames{fileItr});
        [~, fileName, ~] = fileparts(testFileNames{fileItr});
        if max(size(img)) > options.maxImageDim
           img = imresize(img, options.maxImageDim/max(size(img)), 'bilinear'); 
        end
        imwrite(img, [options.processedFolder '/' fileName '.png']);

        %% Form the first level nodes.
        nodes = getNodes(img, options.numberOfFilters, options.gaborFilterThr, ... 
            options.gaborAreaMinResponse, options.gaborFilterSize, currentPath);
        % Keep nodes in the array.
        allNodes((nodeCounter + 1):(nodeCounter + size(nodes,1)), 1:2) = nodes;
        % Assign nodes their image ids.
        imageIds = ones(size(nodes,1), 1)*fileItr;
        allNodes((nodeCounter + 1):(nodeCounter + size(nodes,1)), 3) = ...
                                        mat2cell(imageIds, ones(size(imageIds)));
        % Increment node counter.
        nodeCounter = nodeCounter + size(nodes,1);
    end
    allNodes = allNodes(1:nodeCounter,:);
    graphLevel(size(allNodes,1)) = struct('labelId', [], 'imageId', [], 'position', [], 'children', [], 'parents', [], 'adjInfo', [], 'leafNodes', []);
    nodeImageIds = cell2mat(allNodes(:,3));

    %% Get edges depending on the property to be embedded in the graph.
    [~, edges] = extractEdges(allNodes, [], leafNodeAdjArr, options, 1, options.datasetName, modes);
    %% Fill the basic info in this scene graph level.
    for instanceItr = 1:size(allNodes,1)
       graphLevel(instanceItr).labelId = allNodes{instanceItr,1};
       graphLevel(instanceItr).position = fix(allNodes{instanceItr,2});
       graphLevel(instanceItr).imageId = nodeImageIds(instanceItr);
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
    
    %% Iteratively process each level to parse the object.
    for levelItr = 2:numel(vocabulary)
        graphFileName = [options.testGraphFolder '/' options.datasetName '_' num2str(levelItr-1) '.g'];
        
        %% Visualize the test images with previous layer's subs.
        if options.debug
            visualizeImages( testFileNames, mainGraph, levelItr-1, options, options.datasetName, 'test' );
        end
        
        %% Print the graphs to the input file.
        imageIds = cell2mat(allNodes(:,3));
        numberOfImages = max(imageIds);
        nodeOffset = 0;
        fp = fopen(graphFileName, 'w');
        for imageItr = 1:numberOfImages
            %% Get only nodes and edges belonging to this image and print them.
            imageNodeIdx = find(imageIds==imageItr);
            if numel(edges)>0
                firstNodesOfEdges = edges(:,1);
                imageEdgeIdx = ismember(firstNodesOfEdges, imageNodeIdx);
                imageNodes = allNodes(imageIds==imageItr,:);
                imageEdges = edges(imageEdgeIdx, :);
                imageEdges(:,1:2) = imageEdges(:,1:2) - nodeOffset;
            else
                imageEdges = [];
            end
            fprintf(fp, 'XP\n');
            printGraphToFile(fp, imageNodes(:,1), imageEdges, true);
            nodeOffset = nodeOffset + size(imageNodes,1);
        end
        fclose(fp);
        
        %% Here, we run SUBDUE over the input graph(s) to find pre-defined compositions within the graph.
        % Each pre-defined sub is searched separately.
        newLevel = collectInstances(vocabulary{levelItr}, graphFileName, options);
        
        %% Assign positions, image ids, and leaf nodes. 
        % TODO: It's a shame that we are duplicating code for the following
        % ~100 lines. Fix it!
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
        
        %% Sort newLevel based on the image ids.
        [~, sortedIdx] = sort([newLevel.imageId]);
        newLevel = newLevel(1,sortedIdx);
        
        %% Inhibition! We process the current level to eliminate some of the nodes in the final graph.
        % The rules here are explained in the paper. Basically, each node
        % should introduce a novelty (cover a new area of the image). If it
        % fails to do so, the one which has a lower mdl ranking is
        % discarded. *Natural selection*
%        [newLevel] = applyLocalInhibition(newLevel, options, levelItr);
        
        %% If new level is empty, break.
        if isempty(newLevel)
            break;
        end
        
        %% Create new graph to be fed to for knowledge discovery in the next level.
        numberOfNodes = numel(newLevel);
        allNodes = cell(numel(newLevel), 3);
        for nodeItr = 1:numberOfNodes
           allNodes(nodeItr,:) = {newLevel(nodeItr).labelId, newLevel(nodeItr).position, newLevel(nodeItr).imageId};
        end
        
        %% Create parent relationships.
        mainGraph = mergeIntoGraph(mainGraph, newLevel, levelItr, 1);
        
        %% Extract the edges to form the new graph.
        if numel(modes)<levelItr
           edges=[];
        else
            [~, edges] = extractEdges(allNodes, mainGraph, leafNodeAdjArr, options, levelItr, options.datasetName, modes);
        end
        
        %% If the edges are not empty, fill the edge information in current level.
        if ~isempty(edges)
            for instanceItr = 1:size(allNodes,1)
               nodeEdges = ismember(edges(:,1:2), instanceItr);

               % get non-zero rows of edges
               selfEdges = edges(nodeEdges(:,1) | nodeEdges(:,2),:);
               if numel(selfEdges)>0
                    newLevel(instanceItr).adjInfo = selfEdges;
               end
            end
        end
        
        %% Visualize the test images with previous layer's subs.
        if options.debug && levelItr == numel(vocabulary)
            visualizeImages( testFileNames, mainGraph, levelItr, options, options.datasetName, 'test' );
        end
    end
    filledIdx = find(~(cellfun('isempty', mainGraph)), 1, 'first');
    mainGraph = mainGraph(1:filledIdx);
    save([currentPath '/output/' options.datasetName '_test.mat'], 'mainGraph');
end

