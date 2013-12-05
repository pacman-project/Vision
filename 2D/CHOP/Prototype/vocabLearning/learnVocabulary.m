%> Name: learnVocabulary
%>
%> Description: Given the graph of first-level responses and their
%> relations, this function extracts a hierarchical vocabulary by detecting
%> structural patterns inherent in the input graph. The pattern discovery
%> phase is hierarchic, and continues until a the graph is compressed into
%> a single node.
%>
%> @param allNodes The nodes of the level-1 input graph.
%> @param allEdges The edges of the level-1 input graph.
%> @param graphFileName The input graph's path.
%> @param resultFileName The file which will contain Subdue's output at
%> each level.
%> @param options Program options
%> @param imageIds The image ids of each node.
%> @param currentFolder The path to the workspace folder
%>
%> @retval vocabulary The hierarchic vocabulary learnt from the data.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 26.11.2013
%> Ver 1.1 on 02.12.2013 Completed unlimited number of graph generation.
function [ vocabulary ] = learnVocabulary( allNodes, allEdges, graphFileName, ...
                                                            resultFileName,...
                                                            options, imageIds, currentFolder)
    % Allocate space for vocabulary and hierarchical graph structure
    vocabulary = cell(options.maxLevels,1);
    mainGraph = cell(options.maxLevels,1);
    
    %% Create first vocabulary and graph layers with existing node/edge info
    % Allocate space for current vocabulary level.
    vocabLevel(options.numberOfFilters) = struct('label', [], 'children', [], 'parents', [], 'adjInfo', []);
    % Allocate space for current graph level.
    graphLevel(size(allNodes,1)) = struct('labelId', [], 'imageId', [], 'position', [], 'children', [], 'parents', [], 'adjInfo', []);
    
    % Fill the first layer information
    for subItr = 1:options.numberOfFilters
       vocabLevel(subItr).label = num2str(subItr);
    end
    for instanceItr = 1:size(allNodes,1)
       graphLevel(instanceItr).labelId = allNodes{instanceItr,1};
       graphLevel(instanceItr).position = fix(allNodes{instanceItr,2});
       graphLevel(instanceItr).imageId = imageIds(instanceItr);
       edges = ismember(allEdges(:,1:2), instanceItr);
       
       % get non-zero rows of edges
       selfEdges = allEdges(edges(:,1) | edges(:,2),:);
       if numel(selfEdges)>0
            graphLevel(instanceItr).adjInfo = selfEdges;
       end
    end
    vocabulary(1) = {vocabLevel};
    mainGraph(1) = {graphLevel};
    
    for levelItr = 2:options.maxLevels
        %% Run SUBDUE on the graph for the first time to go from level 1 to level 2.
        runSUBDUE(graphFileName, resultFileName, options, currentFolder);
        
        %% Parse the result file to extract nodes and their relations
        [vocabLevel, graphLevel] = parseResultFile(resultFileName, options);
        
        % If no new subs have been found, finish processing.
        if isempty(vocabLevel)
           break; 
        end
        
        %% Analyze compositions in terms of nodes' spatial distribution
        % Here, we create new substructures based on the spatial
        % distribution of samples.
        
        %% Create the parent relationships between current level and previous level.
        vocabulary = mergeIntoGraph(vocabulary, vocabLevel, levelItr, 0);
        mainGraph = mergeIntoGraph(mainGraph, graphLevel, levelItr, 1);
        
        %% Create new graph to be fed to for knowledge discovery in the next level.
        currentLevel = mainGraph{levelItr};
        numberOfNodes = numel(currentLevel);
        newNodes = cell(numel(currentLevel), 2);
        for nodeItr = 1:numberOfNodes
           newNodes(nodeItr,:) = {currentLevel(nodeItr).labelId, currentLevel(nodeItr).position};
        end
        
        imageIds = [currentLevel(:).imageId];
        imageCount = max(imageIds);
        
        fp = fopen(graphFileName, 'w');
        for fileItr = 1:imageCount
            nodes = newNodes(imageIds==fileItr,:);
            edges = getEdges(nodes, options, levelItr);

            % Print the graph to the file.
            printGraphToFile(fp, nodes(:,1), edges, true);
        end
        fclose(fp);
        
        %% Clear data structures
        clear vocabLevel;
        clear graphLevel;
    end
    
    vocabulary = vocabulary(1:levelItr);
end

%> Name: mergeIntoGraph
%>
%> Description: Given the hierarchical graph and newly formed level, this
%> function embeds new level into the hierarchy by forming parent links
%> between current level and the previous level. If the graph structure has
%> position information, the positions of the children are averaged to find
%> out the position of their super-structure.
%>
%> @param graph Hierarchical graph structure.
%> @param level New level.
%> @param levelItr Level iterator.
%> @param position Position calculations are processed if 1.
%>
%> @retval graph Newly formed hierarchical graph structure.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 02.12.2013
function [graph] = mergeIntoGraph(graph, level, levelItr, position)
    %% Go over children list of each instance in current level
    previousLevel = graph{levelItr-1};
    for newInstItr = 1:numel(level)
        children = level(newInstItr).children;
        for childItr = 1:numel(children)
           if ~ismember(newInstItr, previousLevel(children(childItr)).parents)
                previousLevel(children(childItr)).parents = ...
                    [previousLevel(children(childItr)).parents, newInstItr];
           end
        end
        
        % If necessary, process the position calculation/imageId copying work too.
        if position
            newPosition = [0,0];
            for childItr = 1:numel(children)
               newPosition = newPosition + previousLevel(children(childItr)).position;
            end
            newPosition = fix(newPosition / numel(children));
            level(newInstItr).position = newPosition;
            level(newInstItr).imageId = previousLevel(children(1)).imageId;
        end
    end
    
    %% Assign new levels and move on.
    graph(levelItr-1) = {previousLevel};
    graph(levelItr) = {level};
end