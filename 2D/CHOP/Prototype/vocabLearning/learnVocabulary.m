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
%> @param leafNodeAdjArr The adjacency list of leaf nodes.
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
%> Ver 1.2 on 12.01.2014 Commentary changes, for unified code look.
function [ vocabulary, mainGraph, modes ] = learnVocabulary( allNodes, allEdges, firstModes, leafNodeAdjArr, graphFileName, ...
                                                            resultFileName,...
                                                            options, fileList, datasetName)
                                                        
    %% ========== Step 0: Set initial data structures ==========
    vocabulary = cell(options.maxLevels,1);
    mainGraph = cell(options.maxLevels,1);
    imageIds = cell2mat(allNodes(:,3));
    modes = cell(options.maxLevels,1);
    
    %% ========== Step 1: Create first vocabulary and graph layers with existing node/edge info ==========
    
    %% Step 1.1: Allocate space for current vocabulary and graph levels.
    numberOfFilters = getNumberOfFilters(options);
    vocabLevel(numberOfFilters) = struct('label', [], 'children', [], 'parents', [], 'adjInfo', []);
    graphLevel(size(allNodes,1)) = struct('labelId', [], 'imageId', [], 'position', [], 'children', [], 'parents', [], 'adjInfo', [], 'leafNodes', []);
    
    %% Step 1.2: Fill the first layer information with existing data.
    for subItr = 1:numberOfFilters
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
    
    %% Step 1.3: Prepare intermediate data structures for sequential processing.
    vocabulary(1) = {vocabLevel};
    mainGraph(1) = {graphLevel};
    firstLevel = mainGraph{1};
    modes(1) = {firstModes};
    [pathToGraphFile, fileName, ext] = fileparts(graphFileName);
    [pathToResultFile, rFileName, resultExt] = fileparts(resultFileName);
    previousModes = [];
    
    %% ========== Step 2: Infer new parts by discovering frequent subs in data. ==========
    for levelItr = 2:options.maxLevels
        %% Step 2.1: Run knowledge discovery to learn frequent compositions.
        % TOCHANGE: These two functions should be changed to accommodate
        % other knowledge discovery mechanisms. 
        discoverSubs(graphFileName, resultFileName, options, options.currentFolder, []);
        [vocabLevel, graphLevel] = parseResultFile(resultFileName, options);
       
        % If no new subs have been found, finish processing.
        if isempty(vocabLevel)
            %% Write previous level's appearances to the output folder.
           if options.debug
               if strcmp(options.property, 'mode')
                    visualizeLevel( vocabulary{levelItr-1}, levelItr-1, previousModes, options.currentFolder, options, datasetName);
               end
               visualizeImages( fileList, mainGraph, levelItr-1, options, datasetName, 'train' );
           end
           vocabulary = vocabulary(1:(levelItr-1),:);
           mainGraph = mainGraph(1:(levelItr-1),:);
           modes = modes(1:(levelItr-1),:);
           break; 
        end
       
        %% Step 2.2: Assign realizations R of next graph level (l+1), and fill in their bookkeeping info.
        previousLevel = mainGraph{levelItr-1};
        newImageIds = zeros(numel(graphLevel),1);
        for newNodeItr = 1:numel(graphLevel)
            leafNodes = [];
            position = [0,0];
            graphLevel(newNodeItr).imageId = previousLevel(graphLevel(newNodeItr).children(1)).imageId;
            newImageIds(newNodeItr) = graphLevel(newNodeItr).imageId;
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
        
        % Rearrange graph level so it is sorted by image id.
        [~, sortedIdx] = sort(newImageIds);
        graphLevel = graphLevel(1,sortedIdx);
        
        %% Step 2.3: Inhibition! We process the current level to eliminate some of the nodes in the final graph.
        % The rules here are explained in the paper. Basically, each node
        % should introduce a novelty (cover a new area of the image). If it
        % fails to do so, the one which has a lower mdl ranking is
        % discarded. *Natural selection*
        [graphLevel] = applyLocalInhibition(graphLevel, options, levelItr);
        [remainingComps, ~, IC] = unique([graphLevel.labelId], 'stable');
        
        % Eliminate unused compositions from vocabulary.
        vocabLevel = vocabLevel(1, remainingComps);
        
        % Assign new labels of the remaining realizations.
        for newNodeItr = 1:numel(graphLevel)
            graphLevel(newNodeItr).labelId = IC(newNodeItr);
        end
        %% Step 2.4: Create the parent relationships between current level and previous level.
        vocabulary = mergeIntoGraph(vocabulary, vocabLevel, levelItr, 0);
        mainGraph = mergeIntoGraph(mainGraph, graphLevel, levelItr, 1);
        
        %% Step 2.5: Create object graphs G_(l+1) for the next level, l+1.
        % First, get a list of the new nodes.
        currentLevel = mainGraph{levelItr};
        numberOfNodes = numel(currentLevel);
        newNodes = cell(numel(currentLevel), 3);
        for nodeItr = 1:numberOfNodes
           newNodes(nodeItr,:) = {currentLevel(nodeItr).labelId, currentLevel(nodeItr).position, currentLevel(nodeItr).imageId};
        end
        
        % Extract the edges between new realizations to form the new object graphs.
        [currModes, edges, ~] = extractEdges(newNodes, mainGraph, leafNodeAdjArr, options, levelItr, datasetName, []);
        
        % Learn image ids and number of total images.
        imageIds = [currentLevel(:).imageId]';
        imageCount = max(imageIds);
        
        % Set variables for next level discovery process.
        graphFileName = [pathToGraphFile '/' fileName '_' num2str(levelItr) ext];
        resultFileName = [pathToResultFile '/' rFileName '_' num2str(levelItr) resultExt];
        fp = fopen(graphFileName, 'w');
        nodeOffset = 0;
        
        % Assign prolonging data structures
        currentVocabLevel = vocabulary{levelItr};
        previousModes = modes{levelItr-1};
        
        %% Step 2.6: If no new edges found, kill program.
        if isempty(edges)
           vocabulary = vocabulary(1:(levelItr),:);
           mainGraph = mainGraph(1:(levelItr),:);
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
               
               %% Print out the vocabulary words and object images overlaid with words.
               if strcmp(options.property, 'mode')
                   visualizeLevel( currentVocabLevel, levelItr, previousModes, options.currentFolder, options, datasetName);
               end
               % TODO Fix image id missing problem in the last iteration!.
               visualizeImages( fileList, mainGraph, levelItr, options, datasetName, 'train' );
           end
           break;
        else
            %% Assign edges to their respective nodes.
            for nodeItr = 1:numel(graphLevel)
               nodeEdges = edges(edges(:,1) == nodeItr | edges(:,2) == nodeItr, :);
               graphLevel(nodeItr).adjInfo = nodeEdges;
            end
            mainGraph{levelItr} = graphLevel;
        end
        
        %% Step 2.7: Print the graphs to a file for the next iteration.
        for fileItr = 1:imageCount
            nodes = newNodes(imageIds==fileItr,:);
            imageNodeIdx = find(imageIds==fileItr);
            imageEdgeIdx = ismember(edges(:,1), imageNodeIdx);
            imageEdges = edges(imageEdgeIdx,:);
            imageEdges(:,1:2) = imageEdges(:,1:2) - nodeOffset;
            
            % Print the graph to the file.    
            % Print positive graph indicator
            fprintf(fp, 'XP\n');
            printGraphToFile(fp, nodes(:,1), imageEdges, true);
            nodeOffset = nodeOffset + size(nodes,1);
        end
        fclose(fp);
        
        %% Step 2.8: Write previous level's appearances to the output folder.
        if options.debug
            if strcmp(options.property, 'mode')
                visualizeLevel( vocabulary{levelItr-1}, levelItr-1, previousModes, options.currentFolder, options, datasetName);
            end
            visualizeImages( fileList, mainGraph, levelItr-1, options, datasetName, 'train' );
        end
        
        % Add current level's modes to the main 'modes' array.
        modes(levelItr) = {currModes};
    end
%    vocabulary = vocabulary(1:levelItr,:);
end