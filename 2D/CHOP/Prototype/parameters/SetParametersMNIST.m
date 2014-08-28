%> Name: SetParametersDataset
%>
%> Description: Dataset-specific parameters for CHOP. This part is
%> separated from the main code so that each dataset would have its own
%> file. Relieves us from the need to change the parameters in each
%> different experiment by saving them for future use. Yay.
%>
%> @param datasetName Name of the dataset to work on. 
%> @param options Program options.
%>
%> @retval options Program options.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 26.08.2014
function [ options ] = SetParametersMNIST( datasetName, options )
    %% ========== DATASET - RELATED PARAMETERS ==========
    options.datasetName = datasetName;
    options.learnVocabulary = 1; % If 1, new vocabulary is learned. 
    options.testImages = 1;      % If 1, the test images are processed.
    options.numberOfGaborFilters = 6; % Number of Gabor filters at level 1.
    options.numberOfLHOPFilters = 6; % Number of Gabor filters at level 1.
    
        %% ========== LOW - LEVEL FILTER PARAMETERS ==========
    options.filterType = 'auto'; % If 'gabor': Steerable Gabor filters used 
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
    options.gabor.inhibitionRadius = floor(options.gaborFilterSize/2);
                                        % The inhibition radius basically 
                                        % defines the half of the cube's
                                        % size in which other weaker
                                        % responses than the seed node will
                                        % be surpressed.
    
    options.autoFilterSize = 9;         % Size (one side) of a autodetected 
                                        % filter.
                                        
    options.auto.inhibitionRadius = floor(options.autoFilterSize/2)-1;
    options.autoFilterThr = 0.5;      % Min response threshold for convolved 
                                       % features, taken as the percentage 
                                       % of max response in each image.
    options.autoFilterCount = 100;      % Number of auto-detected filters.
    options.autoFilterPatchCount = 100000; % Number of random patches used 
                                           % to find auto-detected filters.
    
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
                                       % downsampled. DEFAULT 0.5
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
    options.reconstructionType = 'true'; % 'true': Actual reconstruction at each 
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
    options.vis.nodeReconstructionChildren = 1000;
    options.minIndividualReconstructionLevel = 4;   % Minimum image reconstruction 
                                         % level for individual part
                                         % printing. At least 1.
    options.useReceptiveField = 1;       % If 0, regular graph generation 
                                         % takes place. If 1, receptive
                                         % fields are enforced during
                                         % learning.
    options.receptiveFieldSize = options.gaborFilterSize*5; % DEFAULT 5
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
                                           
    options.subdue.maxTime = 3600;           % Max. number of seconds 'self' 
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
end

