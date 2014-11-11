%> Name: SetParametersCommon
%>
%> Description: Common parameters for CHOP. This file is called when the
%> dataset we are working on does not have its own parameter set.
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
function [ options ] = SetParametersCIFAR10( datasetName, options )
    %% ========== DATASET - RELATED PARAMETERS ==========
    options.datasetName = datasetName;
    options.learnVocabulary = 1; % If 1, new vocabulary is learned. 
    options.testImages = 1;      % If 1, the test images are processed.
    options.numberOfGaborFilters = 6; % Number of Gabor filters at level 1.
    
        %% ========== LOW - LEVEL FILTER PARAMETERS ==========
    options.filterType = 'auto'; % If 'gabor': Steerable Gabor filters used 
                                  % as feature detectors.
                                  % If 'auto': Autodetected features.
                                  % Random patches are clustered to obtain
                                  % a number of unsupervised features.
    options.gaborFilterThr = 0.05; % Min response threshold for convolved features, 
                                  % taken as the percentage of max response 
                                  % in each image.
    options.absGaborFilterThr = 0; % Absolute response threshold for low-level 
                                   % responses. ~80 for natural images 
                                   % (depends on many factors though, including 
                                   % size of the filter).
    options.gaborFilterSize = 10;       % Size of a gabor filter. Please note 
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
                                        % defines the half of the square's
                                        % size in which weaker responses other 
                                        % than the seed node will
                                        % be surpressed.
    options.autoFilterSize = 6;         % Size (one side) of a autodetected 
                                        % filter. Assumed to be NxNxD.
    options.auto.inhibitionRadius = floor(options.autoFilterSize/2)-1;
    options.autoFilterThr = 0.05;       % Min response threshold for convolved 
                                       % features, assigned as this percentage 
                                       % of the max response in each image.
    options.autoFilterCount = 100;      % Number of auto-detected filters.
    options.autoFilterPatchCount = 200000; % Number of random patches used 
                                           % to find auto-detected filters.
    options.auto.stride = 2;           % Stride to use when extracting first-
                                       % level features. Only works in
                                       % auto-filter mode, since gabors are
                                       % extracted using conv2, convolution
                                       % implementation of matlab.
    options.auto.deadFeatureStd = 0.00;                                   
%    options.auto.deadFeatureStd = 0.04; % In case of auto-learned features, 
                                       % some dead features may come up.
                                       % The standard deviation check is
                                       % used to eliminate uniform
                                       % features, assigned as this percentage 
                                       % of the max std dev in filters.
    options.distType = 'rank'; % If 'euc': Euclidean distance 
                                                   % (normalized by number
                                                   % of nonzero pixels)
                                                   % will define the
                                                   % distance between two
                                                   % filters. If 'man',
                                                   % manifold distance to
                                                   % be used. If 'rank',
                                                   % manifold distance is
                                                   % used, however, the
                                                   % features are ranked
                                                   % with their distances,
                                                   % and the ranking is the
                                                   % new distance function.
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
                                       % and relations are examined.

    %% ========== CRUCIAL METHOD PARAMETERS (COMPLEXITY, RELATIONS) ==========
    options.noveltyThr = 0.5;           % The novelty threshold used in the 
                                        % inhibition process. At least this 
                                        % percent of a neighboring node's leaf 
                                        % nodes should be new so that it is 
                                        % not inhibited by another higher-
                                        % valued one.
    options.edgeNoveltyThr = 0.75;       % The novelty threshold used in the 
                                        % edge generation. At least this 
                                        % percent of a neighbor node's leaf 
                                        % nodes should be new so that they 
                                        % are linked in the object graph.
    options.edgeQuantize = 20;         % This parameter is used to quantize 
                                        % edges in a edgeQuantize x edgeQuantize 
                                        % window. As the receptive field
                                        % grows, each relation is scaled
                                        % down to this window, and then
                                        % quantized. 
    options.scaling = 0.67;            % Each successive layer is downsampled 
                                       % with a ratio of 1/scaling. Actually,
                                       % the image coordinates of 
                                       % realizations are NOT downsampled, 
                                       % but the edge radius (thus receptive 
                                       % field size) is multiplied by
                                       % 1/scaling at each level, creating
                                       % the same effect.
                                       % DEFAULT 0.5.
    options.edgeType = 'centroid';     % If 'centroid', downsampling is
                                       % applied at each layer, and edges
                                       % link spatially adjacent (within
                                       % its neighborhood) nodes. (No other
                                       % opts at the moment)
    options.reconstructionType = 'leaf'; % 'true': Replacing leaf nodes with 
                                         % average node image in image visualization.
                                         % 'leaf': Detected leaf nodes will
                                         % be marked on the image.
    options.vis.nodeReconstructionChildren = 1000; % Max number of children
                                         % to be used in creation of the image
                                         % for every node in the vocabulary.
    options.vis.instancePerNode = 9;     % Should be square of a natural number.
    if strcmp(options.filterType, 'auto')
        options.receptiveFieldSize = options.autoFilterSize*4; % DEFAULT 5
    else
        options.receptiveFieldSize = options.gaborFilterSize*4;
    end                                  % Size (one side) of the receptive field at
                                         % first level. Please note that in
                                         % each level of the hierarchy, the
                                         % receptive field size grows by 
                                         % 1/scaling.
    options.maxNodeDegree = 6;        % (N) closest N nodes are linked for 
                                       % every node in the object graphs.
    options.maxImageDim = options.receptiveFieldSize*20; %Max dimension of the 
                                       % images the algorithm will work
                                       % with. If one size of a image in
                                       % the dataset is larger than this
                                       % value, it will be rescaled to fit
                                       % in a square of
                                       % maxImageDim x maxImageDim. Aspect ratio
                                       % will be preserved. Set to a large
                                       % value to avoid rescaling.
    options.edgeRadius = floor(options.receptiveFieldSize/2); % The edge radius 
                                       % for two subs to be 
                                       % determined as neighbors. Centroids
                                       % taken into account. This is how a
                                       % receptive field is implemented.
                                       % edgeRadius grows at every level
                                       % with the same ratio as the
                                       % receptive field.
    options.maxLevels = 20;    % The maximum level count for training.
    options.maxInferenceLevels = 20; % The maximum level count for testing.
    
    %% ========== INFERENCE PARAMETERS ==========
    options.fastInference = true; % If set, faster inference (involves 
                                  % inhibition) is performed.
    options.favorParam = 1;      % Between 1:100, if it increases, category 
                                 % nodes with peaks towards a single category 
                                 % will be favoured more, rather than those 
                                 % with relatively uniform distribution.
                                 % Used in determining the category of a node.
    
    %% ========== KNOWLEDGE DISCOVERY PARAMETERS ==========
    options.subdue.evalMetric = 'mdl';     % Evaluation metric for part 
                                           % selection in SUBDUE.
                                           % 'mdl', 'size' or 'freq'. 
                                           % 'mdl': minimum description length,
                                           % 'size': size times frequency.
                                           % 'freq': only takes frequency of
                                           % a node into account.
    options.subdue.isMDLExact = false;     % If true, exact mdl is calculated.
                                           % Otherwise, approximate mdl is
                                           % calculated (faster).
    options.subdue.mdlNodeWeight = 8;      % Weight of a node in DL calculations 
                                           % in MDL-based evaluation
                                           % metric. Cost of a node =
                                           % labelId (int, 4 byte) + pointer to
                                           % edges (int, 4 byte) = 8.
    options.subdue.mdlEdgeWeight = 9;      % Weight of an edge in DL calculations 
                                           % in MDL-based evaluation
                                           % metric. Cost of an edge =
                                           % edgeLabelId (int, 4 byte) + 
                                           % destinationNode (int,4 byte) + 
                                           % isDirected (byte, 1 byte) = 9.
    options.subdue.maxTime = 600;          % Max. number of seconds subdue is
                                            % allowed to run. Typically
                                            % around 100 (secs) for toy data. 
                                            % You can set to higher values
                                            % (e.g. 3600 secs) for large
                                            % datasets.
    options.subdue.threshold = 0.03; % Theshold for elastic part matching. 
                                    % Can be in [0,1]. 
                                    % 0: Strict matching, 
                                    % (value -> 1) Matching criterion 
                                    % gets looser.
                                    % This similarity threshold is used to
                                    % group parts in the same level
                                    % together in order to increase
                                    % generalization ability of detected
                                    % parts.
    options.subdue.minSize = 2; % Minimum number of nodes in a composition.
    options.subdue.maxSize = 3; % Maximum number of nodes in a composition.
    options.subdue.nsubs = 5000;  % Maximum number of nodes allowed in a level.
    options.subdue.beam = 100;   % Beam length in SUBDUE' search mechanism.
    options.subdue.overlap = false;   % If true, overlaps between a substructure's 
                                     % instances are considered in the
                                     % evaluation of the sub. Otherwise,
                                     % unique (in terms of node sets) instances 
                                     % are taken into account (DEFAULT).
                                     % However, all possible instances are
                                     % returned anyway in order to
                                     % introduce redundancy in the final
                                     % object graphs.
     options.subdue.supervised = false; % If true, graph search is performed over
				          % the whole data. If not, individual categories 
			                  % are searched, and the vocabularies are then 
			                  % combined.
end

