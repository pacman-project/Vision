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
    options.maxLevels = 100;    % The maximum level count               
    options.maxLabelLength = 100; % The maximum label name length allowed.
    options.maxNumberOfFeatures = 1000000; % Total maximum number of features.
                                  % The last three are not really parameters, 
                                  % put here to avoid hard-coding.
    options.subdue.instanceIndicator = '  Instance ';
    options.subdue.labelIndicator = 'Label: ';
                                % The instance/label indicators are put
                                % here for compliance with subdue output,
                                % need to be changed if subdue output
                                % format is changed.
    options.subdue.minSize = 2; % Minimum number of nodes in a composition 
    options.subdue.maxSize = 4; % Maximum number of nodes in a composition
    options.subdue.nsubs = 100; % Maximum number of nodes allowed in a level
    options.subdue.diverse = 1; % 1 if diversity is forced, 0 otw
    options.subdue.beam = 10;   % Beam length in SUBDUE
    options.subdue.valuebased = 1; % 1 if value-based queue is used, 0 otw
    options.subdue.overlap = 1; % 1 if overlapping instances allowed, 0 otw
    
    % Learn dataset path relative to this m file
    currentFileName = mfilename('fullpath');
    [currentFilePath, ~, ~] = fileparts(currentFileName);
    datasetFolder = [currentFilePath '/datasets/' datasetName '/'];
    addpath([currentFilePath '/utilities']);
    addpath([currentFilePath '/graphTools']);
    
    % Specify name of the graph files
    graphFileName = [currentFilePath '/graphs/' datasetName '.g'];
    resultFileName = [currentFilePath '/output/' datasetName '.txt'];
    fp = fopen(graphFileName, 'w');

    % Get all images under the dataset
    fileNames = fuf([datasetFolder '*_crop.png'], 1, 'detail');
    
    % Extract graphs of each image and add to the final graph
 %   for fileItr = 1:size(fileNames,1) 
    nodeCounter = 0;
    edgeCounter = 0;
    allNodes = cell(options.maxNumberOfFeatures,2);
    allEdges = zeros(options.maxNumberOfFeatures,2);
 
    for fileItr = 1:20
        img = imread(fileNames{fileItr});
        nodes = getNodes(img, options.numberOfFilters, currentFilePath);
        edges = getEdges(nodes, options.property);
        
        % Keep nodes and edges 
        allNodes((nodeCounter + 1):(nodeCounter + size(nodes,1)), :) = nodes;
        allEdges((edgeCounter + 1):(edgeCounter + size(edges,1)), :) = edges + edgeCounter;
    
        % Print the graph to the file.
        printGraphToFile(fp, cell2mat(nodes(:,1)), edges, true);
        nodeCounter = nodeCounter + size(nodes,1);
        edgeCounter = edgeCounter + size(edges,1);
    end
    
    % Learn the vocabulary in an unsupervised manner from the input graphs.
    [vocabulary] = learnVocabulary(allNodes, allEdges, graphFileName, ...
                                    resultFileName, options);
    
    fclose(fp);
end

