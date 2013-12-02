%> Name: runExperiment
%>
%> Description: The entry function to libHOP code. Calls libHOP with
%> desired parameters on specified dataset. For each image in the dataset,
%> we create a graph out of level 1 gabor filter responses, and form
%> edges. Then, we compress samples of each category with SUBDUE to get 
%>
%> @param datasetName Name of the dataset to work on. 
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 18.11.2013

function [] = runExperiment( datasetName )

    % Initial parameters
    options.numberOfFilters = 6; % Number of Gabor filters at level 1
    options.property = 'co-occurence'; % Geometric property to be examined
                                       % 'co-occurence' or
    options.scaling = 0.95;            % Each successive layer is downsampled 
                                       % with a ratio of 1/scaling. Changes
                                       % formation of edges in upper
                                       % layers, since edge radius
                                       % stays the same while images are 
                                       % downsampled.
    options.edgeRadius = 5;            % The edge radius for two subs to be 
                                       % determined as neighbors. Centroids
                                       % taken into account.
    options.maxLevels = 3;    % The maximum level count               
    options.maxLabelLength = 100; % The maximum label name length allowed.
    options.maxNumberOfFeatures = 1000000; % Total maximum number of features.
                                  % The last three are not really parameters, 
                                  % put here to avoid hard-coding.
    options.subdue.instanceIndicator = sprintf('\n  Instance ');
    options.subdue.labelIndicator = sprintf('Label: ');
    options.subdue.nodeIndicator = sprintf('\nv ');
    options.subdue.edgeIndicator = sprintf('\nu ');
    options.subdue.endLineIndicator = sprintf('\n');
    options.subdue.subPrefix = 'SUB_';
                                % The indicators are put
                                % here for compliance with SUBDUE output,
                                % need to be changed if SUBDUE output
                                % format is changed. They are not
                                % parameters, and should not be changed
                                % unless SUBDUE output format is changed.
    options.subdue.minSize = 2; % Minimum number of nodes in a composition 
    options.subdue.maxSize = 4; % Maximum number of nodes in a composition
    options.subdue.nsubs = 20; % Maximum number of nodes allowed in a level
    options.subdue.diverse = 1; % 1 if diversity is forced, 0 otw
    options.subdue.beam = 4;   % Beam length in SUBDUE
    options.subdue.valuebased = 1; % 1 if value-based queue is used, 0 otw
    options.subdue.overlap = 0; % 1 if overlapping instances allowed, 0 otw
    
    % Learn dataset path relative to this m file
    currentFileName = mfilename('fullpath');
    [currentPath, ~, ~] = fileparts(currentFileName);
    datasetFolder = [currentPath '/datasets/' datasetName '/'];
    addpath([currentPath '/utilities']);
    addpath([currentPath '/graphTools']);
    
    % Specify name of the graph files
    graphFileName = [currentPath '/graphs/' datasetName '.g'];
    resultFileName = [currentPath '/output/' datasetName '.txt'];
    fp = fopen(graphFileName, 'w');

    % Get all images under the dataset
    fileNames = fuf([datasetFolder '*_crop.png'], 1, 'detail');
    
    %% Extract graphs of each image as first level of the hierarchy
    nodeCounter = 0;
    edgeCounter = 0;
    allNodes = cell(options.maxNumberOfFeatures,2);
    allEdges = zeros(options.maxNumberOfFeatures,3);
    imageIds = [];
    
%   for fileItr = 1:size(fileNames,1) 
    for fileItr = 1:20
        img = imread(fileNames{fileItr});
        nodes = getNodes(img, options.numberOfFilters, currentPath);
        imageIds = [imageIds, ones(1, size(nodes,1))*fileItr];
        edges = getEdges(nodes, options, 1);
        
        % Keep nodes and edges 
        allNodes((nodeCounter + 1):(nodeCounter + size(nodes,1)), :) = nodes;
        allEdges((edgeCounter + 1):(edgeCounter + size(edges,1)), 1:2) = edges(:,1:2) + edgeCounter;
        allEdges((edgeCounter + 1):(edgeCounter + size(edges,1)), 3) = edges(:,3);
    
        % Print the graph to the file.
        printGraphToFile(fp, nodes(:,1), edges, true);
        nodeCounter = nodeCounter + size(nodes,1);
        edgeCounter = edgeCounter + size(edges,1);
    end
    allNodes = allNodes(1:nodeCounter,:);
    allEdges = allEdges(1:edgeCounter,:);
    fclose(fp);
    
    %% Learn the vocabulary in an unsupervised manner from the input graphs.
    [vocabulary] = learnVocabulary(allNodes, allEdges, graphFileName, ...
                                    resultFileName, options, imageIds, currentPath);
    save([currentPath '/output/' datasetName '_vb.mat'], 'vocabulary');
end

