%> Name: SetParameters
%>
%> Description: The parameter setting function of CHOP. All program
%> parameters are to be set here for a tidy codebase.
%>
%> @param datasetName Name of the dataset to work on. 
%> @retval options Program options.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 10.01.2014
function [ options ] = SetParameters( datasetName )
    options.debug = 1;           % If debug = 1, additional output will be 
                                 % generated to aid debugging process.
    options.datasetName = datasetName;
    options.learnVocabulary = 1; % If 1, new vocabulary is learned. 
    options.testImages = 0;      % If 1, the test images are processed.
    options.numberOfGaborFilters = 6; % Number of Gabor filters at level 1.
    options.numberOfLHOPFilters = 6; % Number of Gabor filters at level 1.
    options.numberOfAutoFilters = 100; % Number of Gabor filters at level 1.
    options.numberOfTrainingImages = 70; % Number of training images to be 
                                 % used in unsupervised vocabulary learning.
    options.numberOfTestImages = 70; % Number of training images to be 
    options.numberOfVocabImagesPerCategory = 2; % Number of vocabulary images 
                                 % to be used in vocabulary learning.
                                 % used in unsupervised vocabulary learning.   
    options.filterType = 'gabor'; % If 'gabor': Steerable Gabor filters used 
                                  % as feature detectors.
                                  % If 'lhop': Steerable Gabor filters in LHOP 
                                  % are used as feature detectors.
                                  % If 'auto': Autodetected features.
                                  % Random patches are clustered to obtain
                                  % a number of unsupervised features.
    options.gaborFilterThr = 0.2; % Response threshold for convolved features, 
                                  % taken as the percentage of max response 
                                  % in each image.
    options.gaborAreaMinResponse = 0.2; % The threshold to define the minimum response 
                                        % of a filter. Lower-valued responses 
                                        % are inhibited in each response's 
                                        % filter area. Threshold is
                                        % calculated relative to the
                                        % maximum response in that image.
    options.gaborFilterSize = 15;       % Size of a gabor filter. Please note 
                                        % that the size also depends on the 
                                        % filter parameters, so consider them 
                                        % all when you change this!
    options.gabor.sigma = 1;            % Gabor filter parameters
    options.gabor.theta = 0;
    options.gabor.lambda = 1.2;
    options.gabor.psi = 0;
    options.gabor.gamma = 0.4;
    
    options.autoFilterSize = 8;         % Size (one side) of a autodetected 
                                        % filter.
    options.autoFilterCount = 100;      % Number of auto-detected filters.
    options.autoFilterPatchCount = 100000; % Number of random patches used 
                                           % to find auto-detected filters.
    options.noveltyThr = 0.5;           % The novelty threshold used in the 
                                        % inhibition process. At least this 
                                        % percent of a neighbor node's leaf 
                                        % nodes should be new.
    options.edgeNoveltyThr = 0.5;       % The novelty threshold used in the 
                                        % inhibition process. At least this 
                                        % percent of a neighbor node's leaf 
                                        % nodes should be new.
    options.property = 'hist'; % Geometric property to be examined
                                       % 'co-occurence': uniform edges 
                                       % 'mode': cluster relative positions
                                       % 'hist': divide space into x.
    options.scaling = 0.75;            % Each successive layer is downsampled 
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
    options.maxNodeDegreeLevel1 = 5;
    options.maxNodeDegree = 5;         % (N) closest N nodes are considered at
                                       % level 1-l.
    options.maxImageDim = options.gaborFilterSize*25;
    options.maximumModes = 10;          % Maximum number of modes allowed for 
                                       % a node pair.
    options.edgeRadius = options.gaborFilterSize*3; % The edge radius for two subs to be 
                                       % determined as neighbors. Centroids
                                       % taken into account.
    options.maxLevels = 10;    % The maximum level count               
    options.maxLabelLength = 100; % The maximum label name length allowed.
    options.maxNumberOfFeatures = 1000000; % Total maximum number of features.
                                  % The following are not really parameters, 
                                  % put here to avoid hard-coding.
    options.maxNumberOfEdges = 1000000;
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
    options.subdue.maxSize = 4; % Maximum number of nodes in a composition
    options.subdue.nsubs = 4000;  % Maximum number of nodes allowed in a level
    options.subdue.diverse = 1; % 1 if diversity is forced, 0 otw
    options.subdue.beam = 2000;   % Beam length in SUBDUE
    options.subdue.valuebased = 1; % 1 if value-based queue is used, 0 otw
    options.subdue.overlap = 1; % 1 if overlapping instances allowed, 0 otw
    options.subdue.winSep = '\'; % If windows, we replace '/' in command line
                                 % with this.
    
    % Learn dataset path relative to this m file
    currentFileName = mfilename('fullpath');
    [currentPath, ~, ~] = fileparts(currentFileName);
    
    % Set folder parameters.
    options.currentFolder = currentPath;
    options.processedFolder = [currentPath '/output/' datasetName '/original'];
    options.outputFolder = [currentPath '/output/' datasetName];
    options.testOutputFolder = [options.outputFolder '/test'];
    options.preDefinedFolder = [currentPath '/output/' datasetName '/preDefined'];
    options.testGraphFolder = [currentPath '/graphs/' datasetName '/test'];
    options.trainGraphFolder = [currentPath '/graphs/' datasetName '/train'];
end

