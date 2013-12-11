%> Name: runExperiment
%>
%> Description: The entry function to libHOP code. Calls libHOP with
%> desired parameters on specified dataset. For each image in the dataset,
%> we create a graph out of level 1 gabor filter responses, and form
%> edges. Then, we compress samples of each category with SUBDUE to get 
%>
%> @param datasetName Name of the dataset to work on. 
%> @param imageExtension The extension of the files to work on. Examples
%> include '.jpg', '.png', '_crop.png'...
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 18.11.2013
%> Ver 1.1 on 05.12.2013 Various parameter additions, 'mode' changes

function [] = runExperiment( datasetName, imageExtension )
    %% Step 0: Set program parameters and global variables.
    options.debug = 1;           % If debug = 1, additional output will be 
                                 % generated to aid debugging process.
    options.numberOfFilters = 6; % Number of Gabor filters at level 1.
    options.numberOfTrainingImages = 15; % Number of training images to be 
                                 % used in unsupervised vocabulary learning.
    options.numberOfTestImages = 5; % Number of training images to be 
                                 % used in unsupervised vocabulary learning.                             
    options.gaborFilterThr = 0.6; % Response threshold for Gabor filter 
                                 % convolution.
    options.gaborAreaMinResponse = 0.1; % The threshold to define the minimum response 
                                        % of a filter. Lower-valued responses 
                                        % are inhibited in each response's 
                                        % filter area.
    options.noveltyThr = 0.5;           % The novelty threshold used in the 
                                        % inhibition process. At least this 
                                        % percent of a neighbor node's leaf 
                                        % nodes should be new.
    options.gaborFilterSize = 15;       % Size of a gabor filter. Please note 
                                        % that the size also depends on the 
                                        % filter parameters, so consider them 
                                        % all when you change this!
    options.property = 'mode'; % Geometric property to be examined
                                       % 'co-occurence' or
    options.scaling = 0.75;            % Each successive layer is downsampled 
                                       % with a ratio of 1/scaling. Changes
                                       % formation of edges in upper
                                       % layers, since edge radius
                                       % stays the same while images are 
                                       % downsampled.
    options.maxImageDim = options.gaborFilterSize*20;
    options.maximumModes = 8;          % Maximum number of modes allowed for 
                                       % a node pair.
    options.edgeRadius = options.gaborFilterSize; % The edge radius for two subs to be 
                                       % determined as neighbors. Centroids
                                       % taken into account.
    options.maxLevels = 10;    % The maximum level count               
    options.maxLabelLength = 100; % The maximum label name length allowed.
    options.maxNumberOfFeatures = 1000000; % Total maximum number of features.
                                  % The last three are not really parameters, 
                                  % put here to avoid hard-coding.
    options.localInhibitionOverlapThr = 0.5; % The maximum allowed threshold 
                                  % of seed node and adjacent node at each 
                                  % layer. Only this portion of their sub-parts 
                                  % can overlap, otherwise second one is 
                                  % discarded.
    options.subdue.instanceIndicator = sprintf('\n  Instance ');
    options.subdue.labelIndicator = sprintf('Label: ');
    options.subdue.nodeIndicator = sprintf('\nv ');
    options.subdue.edgeIndicator = sprintf('\nu ');
    options.subdue.directedEdgeIndicator = sprintf('\nd ');
    options.subdue.endLineIndicator = sprintf('\n');
    options.subdue.subPrefix = 'SUB_';
                                % The indicators are put
                                % here for compliance with SUBDUE output,
                                % need to be changed if SUBDUE output
                                % format is changed. They are not
                                % parameters, and should not be changed
                                % unless SUBDUE output format is changed.
    options.subdue.minSize = 2; % Minimum number of nodes in a composition 
    options.subdue.maxSize = 5; % Maximum number of nodes in a composition
    options.subdue.nsubs = 100;  % Maximum number of nodes allowed in a level
    options.subdue.diverse = 1; % 1 if diversity is forced, 0 otw
    options.subdue.beam = 50;   % Beam length in SUBDUE
    options.subdue.valuebased = 1; % 1 if value-based queue is used, 0 otw
    options.subdue.overlap = 0; % 1 if overlapping instances allowed, 0 otw
                                % Right now, there is a bug with SUBDUE
                                % code, so a value of 0 already allows
                                % overlap and works better. Not advised to
                                % change it at the moment.
    
    % Learn dataset path relative to this m file
    currentFileName = mfilename('fullpath');
    [currentPath, ~, ~] = fileparts(currentFileName);
    % Set folder variables.
    options.currentFolder = currentPath;
    datasetFolder = [currentPath '/datasets/' datasetName '/'];
    options.processedFolder = [currentPath '/output/' datasetName '/original'];
    options.outputFolder = [currentPath '/output/' datasetName];
    
    % Add relevant folders to path.
    addpath([currentPath '/utilities']);
    addpath([currentPath '/graphTools']);
    addpath([currentPath '/vocabLearning']);
    
    % Specify name of the graph files
    graphFileName = [currentPath '/graphs/' datasetName '.g'];
    resultFileName = [options.outputFolder '/' datasetName '.txt'];
    fp = fopen(graphFileName, 'w');
    
    if ~exist(options.processedFolder,'dir')
       mkdir(options.processedFolder); 
    end

    % Get all images under the dataset
    fileNames = fuf([datasetFolder '*', imageExtension], 1, 'detail');
    trainingFileNames = fileNames(1:options.numberOfTrainingImages,:);
    
    %% Step 1: Extract nodes of each image for the first level of the hierarchy.
    nodeCounter = 0;
    allNodes = cell(options.maxNumberOfFeatures,3);
    
   for fileItr = 1:size(trainingFileNames,1) 
        %% First, downsample the image if it is too big.
        img = imread(trainingFileNames{fileItr});
        [~, fileName, ~] = fileparts(trainingFileNames{fileItr});
        if max(size(img)) > options.maxImageDim
           img = imresize(img, options.maxImageDim/max(size(img))); 
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
    
    %% Step 2: Get edges depending on the property to be embedded in the graph.
    [modes, edges] = extractEdges(allNodes, options, 1, datasetName);
    
    %% Step 3: Print the graphs to the input file.
    imageIds = cell2mat(allNodes(:,3));
    numberOfImages = max(imageIds);
    nodeOffset = 0;
    for imageItr = 1:numberOfImages
        %% Get only nodes and edges belonging to this image and print them.
        imageNodeIdx = find(imageIds==imageItr);
        firstNodesOfEdges = edges(:,1);
        imageEdgeIdx = ismember(firstNodesOfEdges, imageNodeIdx);
        imageNodes = allNodes(imageIds==imageItr,:);
        imageEdges = edges(imageEdgeIdx, :);
        imageEdges(:,1:2) = imageEdges(:,1:2) - nodeOffset;
        printGraphToFile(fp, imageNodes(:,1), imageEdges, true);
        nodeOffset = nodeOffset + size(imageNodes,1);
    end
    fclose(fp);
    
    %% Step 4: Learn the vocabulary in an unsupervised manner from the input graphs.
    [vocabulary, mainGraph, modes] = learnVocabulary(allNodes, edges, modes, graphFileName, ...
                                    resultFileName, options, trainingFileNames, datasetName);
    save([currentPath '/output/' datasetName '_vb.mat'], 'vocabulary', 'mainGraph', 'modes');
    
    %% Step 5: For each test image, run inference code. 
    % Test images are obtained from the rest of images in dataset.
%    testFileNames = fileNames((options.numberOfTrainingImages+1):(options.numberOfTrainingImages+options.numberOfTestImages),:);
%    preparePreDefinedFiles(vocabulary);
    
end

