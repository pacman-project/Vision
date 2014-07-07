%> Name: SetParameters
%>
%> Description: The parameter setting function of CHOP. All program
%> parameters are to be set here for a tidy codebase. In addition, some
%> initializations and folder generations also run here.
%>
%> @param datasetName Name of the dataset to work on. 
%>
%> @retval options Program options.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 10.01.2014
function [ options ] = SetParameters( datasetName, isTraining )
    options.isTraining = isTraining;
    %% ========== DEBUG PARAMETER ==========
    options.debug = 1;           % If debug = 1, additional output will be 
                                 % generated to aid debugging process.
                                 
    %% ========== PARALLEL PROCESSING PARAMETERS ==========
    options.parallelProcessing = true;
    options.numberOfThreads = feature('NumCores') * 2 - 1;
    %% ========== DATASET - RELATED PARAMETERS ==========
    options.datasetName = datasetName;
    options.learnVocabulary = 1; % If 1, new vocabulary is learned. 
    options.testImages = 1;      % If 1, the test images are processed.
    options.numberOfGaborFilters = 6; % Number of Gabor filters at level 1.
    options.numberOfLHOPFilters = 6; % Number of Gabor filters at level 1.
    
    % Rest of the parameters in this section are obsolete, please ignore
    % them.
    options.numberOfAutoFilters = 100; % Number of Gabor filters at level 1.  
        %% ========== LOW - LEVEL FILTER PARAMETERS ==========
    options.filterType = 'gabor'; % If 'gabor': Steerable Gabor filters used 
                                  % as feature detectors.
                                  % If 'lhop': Steerable Gabor filters in LHOP 
                                  % are used as feature detectors.
                                  % If 'auto': Autodetected features.
                                  % Random patches are clustered to obtain
                                  % a number of unsupervised features.
    options.gaborFilterThr = 0.075; % Min response threshold for convolved features, 
                                  % taken as the percentage of max response 
                                  % in each image.
    options.absGaborFilterThr = 0; % Absolute response threshold for low-level 
                                   % responses. ~80 for natural images.
    options.gaborFilterSize = 9;       % Size of a gabor filter. Please note 
                                        % that the size also depends on the 
                                        % filter parameters, so consider them 
                                        % all when you change this!
    options.gabor.sigma = 1;            % Gabor filter parameters
    options.gabor.theta = 0;
    options.gabor.lambda = 1;
    options.gabor.psi = 0;
    options.gabor.gamma = 0.25;
    options.gabor.inhibitionRadius = floor(options.gaborFilterSize/2)+2;
                                        % The inhibition radius basically 
                                        % defines the half of the cube's
                                        % size in which other weaker
                                        % responses than the seed node will
                                        % be surpressed.
    
    options.autoFilterSize = 9;         % Size (one side) of a autodetected 
                                        % filter.
    options.autoFilterCount = 64;      % Number of auto-detected filters.
    options.autoFilterPatchCount = 50000; % Number of random patches used 
                                           % to find auto-detected filters.
    options.autoFilterVisX = 8;        % Visualization parameters for 
                                        % auto-generated filters. The
                                        % filters are visualized in a
                                        % grid-like image specified with
                                        % these parameters.
    options.autoFilterVisY = 8;
    
    %% ========== GT Parameters ==========
    options.useGT = true;              % If true, gt info is used. 
    options.gtType = 'contour';        % 'contour' type gt: nodes lying under
                                       % the gt contour is examined (within
                                       % a neighborhood defined by
                                       % contourGTNeighborhood). 
                                       % 'bbox' type gt: nodes in the
                                       % gt bounding box are examined.
    options.contourGTNeighborhood = 8;% width of the band along the contour 
                                       % (half width, actual width is
                                       % double this value)  in which nodes
                                       % are examined.

    %% ========== CRUCIAL METHOD PARAMETERS (COMPLEXITY, RELATIONS) ==========
    options.noveltyThr = 0.5;           % The novelty threshold used in the 
                                        % inhibition process. At least this 
                                        % percent of a neighboring node's leaf 
                                        % nodes should be new so that it is 
                                        % not inhibited by center.
    options.edgeNoveltyThr = 0.75;       % The novelty threshold used in the 
                                        % edge generation. At least this 
                                        % percent of a neighbor node's leaf 
                                        % nodes should be new so that they 
                                        % are linked in the object graph.
    options.property = 'mode'; % Geometric property to be examined
                                       % 'co-occurence': uniform edges 
                                       % 'mode': clusters of relative positions
                                       % 'hist': divide space into 8 
                                       % pre-defined regions.
    options.mode.maxSamplesPerMode = 200; % In mode calculation between node1
                                          % and node2, not all samples are
                                          % considered. Randomly chosen
                                          % samples are used, defined with
                                          % this number.
    options.mode.minSamplesPerMode = 4;   % The minimum number of samples to 
                                          % be assigned to each mode. If
                                          % there are not enough samples,
                                          % number of modes for that
                                          % specific part pair is reduced
                                          % automatically to match this
                                          % number, if possible.
    options.scaling = 0.5;            % Each successive layer is downsampled 
                                       % with a ratio of 1/scaling. Changes
                                       % formation of edges in upper
                                       % layers, since edge radius
                                       % stays the same while images are 
                                       % downsampled.
    options.edgeType = 'centroid';      % If 'contour', the nodes in upper layers 
                                       % are linked if their leaf nodes are 
                                       % neighbors in the first layer.If
                                       % 'centroid', downsampling is
                                       % applied at each layer, and edges
                                       % link spatially adjacent (within
                                       % its neighborhood) nodes.
                                       % UPDATE: 'contour' TYPE EDGE LACKS
                                       % CORRECT IMPLEMENTATION, USE
                                       % 'centroid' INSTEAD.
    options.highLevelProperty = 'none'; % If high-level edges between compositions
                                       % in different objects' graphs are
                                       % considered, they are specified
                                       % here. 'multiview' and 'category'
                                       % are considered to be implemented,
                                       % among others.
    options.reconstructionType = 'leaf'; % 'true': Actual reconstruction at each 
                                         % level with compositions' masks on 
                                         % that specific level
                                         % 'leaf': Detected leaf nodes will
                                         % be marked on the image.
    options.imageReconstructionType = 'all'; % If 'individual' type 
                                         % image reconstruction is used,
                                         % each realization is written to a
                                         % different image, along with its
                                         % normalized mdl score. If 'all'
                                         % type image reconstruction is
                                         % used, all realizations are
                                         % written in a single image.
    options.minIndividualReconstructionLevel = 4;   % Minimum image reconstruction 
                                         % level for individual part
                                         % printing. At least 1.
    options.useReceptiveField = 1;       % If 0, regular graph generation 
                                         % takes place. If 1, receptive
                                         % fields are enforced during
                                         % learning.
    options.receptiveFieldSize = options.gaborFilterSize*5;
                                         % Size (one side) of the receptive field at
                                         % each level. Please note that in
                                         % each level of the hierarchy, the
                                         % coordinates are downsampled, so our
                                         % receptive field indeed grows.
    options.maxNodeDegreeLevel1 = 7;
    options.maxNodeDegree = 6;         % (N) closest N nodes are considered at
                                       % level 1-l, to link nodes via edges.
                                       % UPDATE: If receptive fields are
                                       % used, no max degree is applied.
%    options.maxImageDim = options.gaborFilterSize*40;
    options.maxImageDim = options.gaborFilterSize*1000;
    options.maximumModes = 50;          % Maximum number of modes allowed for 
                                       % a node pair.
    options.edgeRadius = floor(options.receptiveFieldSize/2); % The edge radius for two subs to be 
                                       % determined as neighbors. Centroids
                                       % taken into account.
    
    options.maxLevels = 20;    % The maximum level count for training.
    options.maxInferenceLevels = 20; % The maximum level count for testing.
                                    % Please write 1 off.
    options.maxLabelLength = 100; % The maximum label name length allowed.
    
    %% ========== INFERENCE PARAMETERS ==========
    options.fastInference = true;
    
    %% ========== KNOWLEDGE DISCOVERY PARAMETERS ==========
                                           % The following metric is valid
                                           % only in 'self' implementation.
    options.subdue.evalMetric = 'mdl';     % 'mdl' or 'size'. 'mdl' takes 
                                           % the relations between
                                           % receptive fields into account,
                                           % while 'size' based metric
                                           % treats each receptive field as
                                           % separate graphs, and evaluates
                                           % subs based on (size x
                                           % frequency).
                                           
    options.subdue.maxTime = 1200;           % Max. number of seconds 'self' 
                                            % type implemented subdue is
                                            % run over data. Typically
                                            % around 100 (secs).
    options.subdue.maxInferenceTime = 0.5;  % Max. number of seconds 'self' 
                                            % type subdue is run over data,
                                            % given a pre-defined
                                            % vocabulary. Typically around
                                            % 0.5 (secs).
                                % The indicators are put
                                % here for compliance with SUBDUE output,
                                % need to be changed if SUBDUE output
                                % format is changed. They are not
                                % parameters, and should not be changed
                                % unless SUBDUE output format is changed.
    options.subdue.threshold = 0.05; % Theshold for elasticity-based matching 
                                    % in SUBDUE. Can be in [0,1]. 0: Strict
                                    % matching, (-> 1) Matching gets looser.
    options.subdue.minSize = 2; % Minimum number of nodes in a composition 
    options.subdue.maxSize = 3; % Maximum number of nodes in a composition
    options.subdue.nsubs = 10000;  % Maximum number of nodes allowed in a level
    options.subdue.beam = 200;   % Beam length in SUBDUE
    options.subdue.winSep = '\'; % If windows, we replace '/' in command line
                                 % with this.
    
    %% ========== FOLDER STRUCTURE INITIALIZATION & GENERATION ==========
    % Learn dataset path relative to this m file
    currentFileName = mfilename('fullpath');
    [currentPath, ~, ~] = fileparts(currentFileName);
    
    % Set folder parameters.
    options.currentFolder = currentPath;
    options.debugFolder = [currentPath '/debug/' datasetName];
    options.processedFolder = [currentPath '/output/' datasetName '/original'];
    options.processedGTFolder = [currentPath '/output/' datasetName '/gt'];
    options.outputFolder = [currentPath '/output/' datasetName];
    options.testOutputFolder = [options.outputFolder '/test'];
    options.testInferenceFolder = [options.outputFolder '/test/inference'];
    options.preDefinedFolder = [currentPath '/output/' datasetName '/preDefined'];
    options.testGraphFolder = [currentPath '/graphs/' datasetName '/test'];
    options.trainGraphFolder = [currentPath '/graphs/' datasetName '/train'];
    options.smoothedFolder = [currentPath '/output/' datasetName '/smoothed'];
    
    %% ========== PATH FOLDER ADDITION ==========
    w = warning('off', 'all');
    addpath(genpath([options.currentFolder '/utilities']));
    addpath(genpath([options.currentFolder '/demo']));
    addpath(genpath([options.currentFolder '/graphTools']));
    addpath(genpath([options.currentFolder '/vocabLearning']));
    addpath(genpath([options.currentFolder '/inference']));
    addpath(genpath([options.currentFolder '/categorization']));
    warning(w);
    
    %% ========== INTERNAL DATA STRUCTURES ==========
    % Internal data structure for a vocabulary level.
    options.vocabNode = VocabNode;
    % Internal data structure for an graph representing a set of object 
    % graphs in a given level.                        
    options.graphNode = GraphNode;     % 1: positive, 0: negative node.
    
    %% ========== LOW - LEVEL FILTER GENERATION ==========
    filters = createFilters(options);
    options.filters = filters;
    options.numberOfFilters = numel(filters);
    
    if strcmp(options.property, 'co-occurence') && strcmp(options.reconstructionType, 'true')
        options.reconstructionType = 'leaf';
        display('"co-occurence" property and "true" reconstruction is incompatible. Switching to "leaf" type reconstruction.');
    end
end

