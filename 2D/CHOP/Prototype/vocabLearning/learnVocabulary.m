%> Name: learnVocabulary
%>
%> Description: Given the graph of first-level responses and their
%> relations, this function extracts a hierarchical vocabulary by detecting
%> structural patterns inherent in the input graph. The pattern discovery
%> phase is hierarchic, and continues until a the graph is compressed into
%> a single node.
%>
%> @param allNodes The nodes of the level 1 input graph. Each node consists
%> of a label id, position and image id.
%> @param allEdges The edges of the level 1 input graph.
%> @param firstModes The pairwise modes in the leve1 input graph.
%> @param graphFileName The input graph's path.
%> @param resultFileName The file which will contain Subdue's output at
%> each level.
%> @param options Program options
%> @param fileList input image name list.
%> @param datasetName The name of the dataset.
%>
%> @retval vocabulary The hierarchic vocabulary learnt from the data.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 26.11.2013
%> Ver 1.1 on 02.12.2013 Completed unlimited number of graph generation.
function [ vocabulary, mainGraph, modes ] = learnVocabulary( allNodes, allEdges, firstModes, graphFileName, ...
                                                            resultFileName,...
                                                            options, fileList, datasetName)
    %% Allocate space for vocabulary and hierarchical graph structure
    vocabulary = cell(options.maxLevels,1);
    mainGraph = cell(options.maxLevels,1);
    imageIds = cell2mat(allNodes(:,3));
    modes = cell(options.maxLevels,1);
    %% Create first vocabulary and graph layers with existing node/edge info
    % Allocate space for current vocabulary level.
    vocabLevel(options.numberOfFilters) = struct('label', [], 'children', [], 'parents', [], 'adjInfo', []);
    % Allocate space for current graph level.
    graphLevel(size(allNodes,1)) = struct('labelId', [], 'imageId', [], 'position', [], 'children', [], 'parents', [], 'adjInfo', [], 'leafNodes', []);
    
    % Fill the first layer information
    for subItr = 1:options.numberOfFilters
       vocabLevel(subItr).label = num2str(subItr);
    end
    for instanceItr = 1:size(allNodes,1)
       graphLevel(instanceItr).labelId = allNodes{instanceItr,1};
       graphLevel(instanceItr).position = fix(allNodes{instanceItr,2});
       graphLevel(instanceItr).imageId = imageIds(instanceItr);
       graphLevel(instanceItr).leafNodes = instanceItr;
       edges = ismember(allEdges(:,1:2), instanceItr);
       
       % get non-zero rows of edges
       selfEdges = allEdges(edges(:,1) | edges(:,2),:);
       if numel(selfEdges)>0
            graphLevel(instanceItr).adjInfo = selfEdges;
       end
    end
    
    %% Prepare intermediate data structures for sequential processing.
    vocabulary(1) = {vocabLevel};
    mainGraph(1) = {graphLevel};
    firstLevel = mainGraph{1};
    modes(1) = {firstModes};
    [pathToGraphFile, fileName, ext] = fileparts(graphFileName);
    [pathToResultFile, rFileName, resultExt] = fileparts(resultFileName);
    previousModes = [];
    
    for levelItr = 2:options.maxLevels
        %% Run SUBDUE on the graph for the first time to go from level 1 to level 2.
        runSUBDUE(graphFileName, resultFileName, options, options.currentFolder);
        
        %% Parse the result file to extract nodes and their relations
        [vocabLevel, graphLevel] = parseResultFile(resultFileName, options);
        
        % If no new subs have been found, finish processing.
        if isempty(vocabLevel)
           vocabulary = vocabulary(1:(levelItr-1),:);
           mainGraph = mainGraph(1:(levelItr-1),:);
           modes = modes(1:(levelItr-1),:);
           break; 
        end
        
        % Assign positions, image ids and leaf nodes.
        previousLevel = mainGraph{levelItr-1};
        for newNodeItr = 1:numel(graphLevel)
            leafNodes = [];
            position = [0,0];
            graphLevel(newNodeItr).imageId = previousLevel(graphLevel(newNodeItr).children(1)).imageId;
            for childItr = 1:numel(graphLevel(newNodeItr).children)
               leafNodes = [leafNodes, previousLevel(graphLevel(newNodeItr).children(childItr)).leafNodes]; 
            end
            leafNodes = unique(leafNodes);
            for leafNodeItr = 1:numel(leafNodes)
               position = position + firstLevel(leafNodes(leafNodeItr)).position;
            end
            graphLevel(newNodeItr).leafNodes = leafNodes;
            graphLevel(newNodeItr).position = round(position/numel(leafNodes));
        end
        
        %% Inhibition! We process the current level to eliminate some of the nodes in the final graph.
        % The rules here are explained in the paper. Basically, each node
        % should introduce a novelty (cover a new area of the image). If it
        % fails to do so, the one which has a lower mdl ranking is
        % discarded. *Natural selection*
        [graphLevel] = applyLocalInhibition(graphLevel, options, levelItr);
        [remainingComps, ~, newLabelIds] = unique([graphLevel.labelId], 'stable');
        
        %% Eliminate unused compositions from vocabulary.
        vocabLevel = vocabLevel(1, remainingComps);
        for newAssgnItr = 1:numel(graphLevel)
            graphLevel(newAssgnItr).labelId = newLabelIds(newAssgnItr);
        end
        
        %% Create the parent relationships between current level and previous level.
        vocabulary = mergeIntoGraph(vocabulary, vocabLevel, levelItr, 0);
        mainGraph = mergeIntoGraph(mainGraph, graphLevel, levelItr, 1);
        
        %% Write previous level's appearances to the output folder.
        if options.debug
           visualizeLevel( vocabulary{levelItr-1}, levelItr-1, previousModes, options.currentFolder, options, datasetName);
           visualizeImages( fileList, mainGraph, levelItr-1, options, datasetName, 'train' );
        end
        
        %% Create new graph to be fed to for knowledge discovery in the next level.
        currentLevel = mainGraph{levelItr};
        numberOfNodes = numel(currentLevel);
        newNodes = cell(numel(currentLevel), 3);
        for nodeItr = 1:numberOfNodes
           newNodes(nodeItr,:) = {currentLevel(nodeItr).labelId, currentLevel(nodeItr).position, currentLevel(nodeItr).imageId};
        end
        
        %% Extract the edges to form the new graph.
        [currModes, edges] = extractEdges(newNodes, options, levelItr, datasetName);
        
        % Learn image ids and number of total images.
        imageIds = [currentLevel(:).imageId]';
        imageCount = max(imageIds);
        
        % Set variables for next level discovery process.
        graphFileName = [pathToGraphFile '/' fileName '_' num2str(levelItr) ext];
        resultFileName = [pathToResultFile '/' rFileName '_' num2str(levelItr) resultExt];
        fp = fopen(graphFileName, 'w');
        nodeOffset = 0;
        
        %% Assign prolonging data structures
        currentLevel = vocabulary{levelItr};
        previousModes = modes{levelItr-1};
        
        %% If no new edges found, kill program.
        if isempty(edges)
           vocabulary = vocabulary(1:(levelItr-1),:);
           mainGraph = mainGraph(1:(levelItr-1),:);
           modes = modes(1:(levelItr-1),:);
           % Print new level's words before exiting from loop
           if options.debug
               
                % Assign positions, image ids and leaf nodes.
                previousLevel = mainGraph{levelItr-1};
                for newNodeItr = 1:numel(graphLevel)
                    leafNodes = [];
                    position = [0,0];
                    graphLevel(newNodeItr).imageId = previousLevel(graphLevel(newNodeItr).children(1)).imageId;
                    for childItr = 1:numel(graphLevel(newNodeItr).children)
                       leafNodes = [leafNodes, previousLevel(graphLevel(newNodeItr).children(childItr)).leafNodes]; 
                    end
                    leafNodes = unique(leafNodes);
                    for leafNodeItr = 1:numel(leafNodes)
                       position = position + firstLevel(leafNodes(leafNodeItr)).position;
                    end
                    graphLevel(newNodeItr).leafNodes = leafNodes;
                    graphLevel(newNodeItr).position = round(position/numel(leafNodes));
               end
               
               
               vocabulary = mergeIntoGraph(vocabulary, vocabLevel, levelItr, 0);
               mainGraph = mergeIntoGraph(mainGraph, graphLevel, levelItr, 1);
               
               visualizeLevel( currentLevel, levelItr, previousModes, options.currentFolder, options, datasetName);
               % TODO Fix image id missing problem in the last iteration!.
               visualizeImages( fileList, mainGraph, levelItr, options, datasetName, 'train' );
           end
           break; 
        end
        
        %% Print the graphs to the file for new level discovery.
        for fileItr = 1:imageCount
            nodes = newNodes(imageIds==fileItr,:);
            imageNodeIdx = find(imageIds==fileItr);
            imageEdgeIdx = ismember(edges(:,1), imageNodeIdx);
            imageEdges = edges(imageEdgeIdx,:);
            imageEdges(:,1:2) = imageEdges(:,1:2) - nodeOffset;
            
            % Print the graph to the file.
            printGraphToFile(fp, nodes(:,1), imageEdges, true);
            nodeOffset = nodeOffset + size(nodes,1);
        end
        fclose(fp);
        
        % Add current level's modes to the main 'modes' array.
        modes(levelItr) = {currModes};
        
        %% Clear data structures
 %       clear vocabLevel;
 %       clear graphLevel;
 %       clear modes;
    end
    vocabulary = vocabulary(1:(levelItr-1),:);
end