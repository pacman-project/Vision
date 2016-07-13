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
function [ options ] = SetParametersCaltech101_filter( datasetName, options )
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
    options.gaborFilterThr = 0.1; % Min response threshold for convolved features, 
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
    options.gabor.inhibitionRadius = floor(options.gaborFilterSize/2)-1;
                                        % The inhibition radius basically 
                                        % defines the half of the square's
                                        % size in which weaker responses other 
                                        % than the seed node will
                                        % be surpressed.
    options.autoFilterSize = 11;         % Size (one side) of a autodetected 
                                        % filter. Assumed to be NxNxD.
    options.auto.inhibitionRadius = floor(options.autoFilterSize/2)-1;
    options.autoFilterThr = 0.1;       % Min response threshold for convolved 
                                       % features, assigned as this percentage 
                                       % of the max response in each image.
    options.autoFilterCount = 100;      % Number of auto-detected filters.
    options.autoFilterPatchCount = 100000; % Number of random patches used 
                                           % to find auto-detected filters.
    options.auto.stride = 4;           % Stride to use when extracting first-
                                       % level features. Only works in
                                       % auto-filter mode, since gabors are
                                       % extracted using conv2, convolution
                                       % implementation of matlab.                                 
    options.auto.deadFeatureStd = 0.05; % In case of auto-learned features, 
                                       % some dead features may come up.
                                       % The standard deviation check is
                                       % used to eliminate uniform
                                       % features, assigned as this percentage 
                                       % of the max std dev in filters.
    options.distType = 'euc'; % If 'euc': Euclidean distance 
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
    options.noveltyThr = 0.4;           % The novelty threshold used in the 
                                        % inhibition process. At least this 
                                        % percent of a neighboring node's leaf 
                                        % nodes should be new so that it is 
                                        % not inhibited by another higher-
                                        % valued one.
    options.edgeNoveltyThr = 0.8;       % The novelty threshold used in the 
                                        % edge generation. At least this 
                                        % percent of a neighbor node's leaf 
                                        % nodes should be new so that they 
                                        % are linked in the object graph.
    options.edgeQuantize = 7;         % This parameter is used to quantize 
                                        % edges in a edgeQuantize x edgeQuantize 
                                        % window. As the receptive field
                                        % grows, each relation is scaled
                                        % down to this window, and then
                                        % quantized. 
    options.scaling = 0.7;            % Each successive layer is downsampled 
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
    options.minContinuityCoverage = 0.9; % If data coverage drops below this,
                                         % we switch to 'centroid' nodes.
    options.reconstructionType = 'leaf'; % 'true': Replacing leaf nodes with 
                                         % average node image in image visualization.
                                         % 'leaf': Detected leaf nodes will
                                         % be marked on the image.
    options.vis.nodeReconstructionChildren = 1000; % Max number of children
                                         % to be used in creation of the image
                                         % for every node in the vocabulary.
    options.vis.instancePerNode = 9;     % Should be square of a natural number.
    options.vis.visualizedNodes = 100; % Number of vocabulary nodes to be visualized.
    if strcmp(options.filterType, 'auto')
        options.receptiveFieldSize = options.autoFilterSize*2.5; % DEFAULT 5
    else
        options.receptiveFieldSize = options.gaborFilterSize*2.5;
    end                                  % Size (one side) of the receptive field at
                                         % first level. Please note that in
                                         % each level of the hierarchy, the
                                         % receptive field size grows by 
                                         % 1/scaling.
    options.maxNodeDegree = 8;        % (N) closest N nodes are linked for 
                                       % every node in the object graphs.
    options.maxImageDim = options.receptiveFieldSize*8; %Max dimension of the 
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
    options.maxLevels = 10;    % The maximum level count for training.
    options.maxInferenceLevels = 10; % The maximum level count for testing.
    
    %% ========== INFERENCE PARAMETERS ==========
    options.fastInference = true; % If set, faster inference (involves 
                                  % inhibition) is performed.
    options.favorParam = 1;      % Between 1:100, if it increases, category 
                                 % nodes with peaks towards a single category 
                                 % will be favoured more, rather than those 
                                 % with relatively uniform distribution.
                                 % Used in determining the category of a node.
                                 
    %% ========== RECONSTRUCTION PARAMETERS ==========
    options.reconstruction.stoppingCoverage = 1; % Between [0.00, 1.00].
                                           % The default value is 0.99.
                                           % When the training data
                                           % coverage is reached to this
                                           % percent, reconstructive part 
                                           % selection stops.
    options.reconstruction.numberOfReconstructiveSubs = 300; % The maximum 
                                           % number of reconstructive parts
                                           % that can be selected.

    %% ========== GRAPH MATCHING PARAMETERS ==========
    options.nodeSimilarityAllowed = false; % If true, node similarities are 
                                           % considered in graph matching.
                                           % If not, identicality in labels
                                           % represents zero-cost matching,
                                           % while every other kind of node
                                           % correspondance yields a cost
                                           % of 1 (max value). 
    options.edgeSimilarityAllowed = true;  % If true, edge similarities are 
                                           % considered in graph matching.
                                           % If not, identicality in labels
                                           % represents zero-cost matching,
                                           % while every other kind of edge
                                           % transformation yields a cost
                                           % of 1 (max value). 

    %% ========== KNOWLEDGE DISCOVERY PARAMETERS ==========
    options.subdue.evalMetric = 'size';     % Evaluation metric for part 
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
    options.subdue.threshold = 0.05; % Theshold for elastic part matching. 
                                    % Can be in [0,1]. 
                                    % 0: Strict matching, 
                                    % (value -> 1) Matching criterion 
                                    % gets looser.
                                    % This similarity threshold is used to
                                    % group parts in the same level
                                    % together in order to increase
                                    % generalization ability of detected
                                    % parts.
                                    % Ignored if reconstruction flag is
                                    % true, since an optimal threshold is
                                    % searched within the limits specified
                                    % by minThreshold and maxThreshold.
    % The following min/max threshold values limit the area in which an
    % optimal elasticity threshold is going to be searched. 
    options.subdue.minThreshold = 0.01; % Minimum threshold for elastic matching.
    options.subdue.maxThreshold = 0.05; % Max threshold for elastic part matching. 
    options.subdue.thresholdSearchMaxDepth = 4; % The depth of binary search 
                                % when looking for an optimal threshold.
    options.subdue.minSize = 1; % Minimum number of nodes in a composition.
    options.subdue.maxSize = 3; % Maximum number of nodes in a composition.
    options.subdue.nsubs = 50000;  % Maximum number of nodes allowed in a level.
    options.subdue.beam = 100;   % Beam length in SUBDUE' search mechanism.
    options.subdue.overlap = false;   % If true, overlaps between a substructure's 
                                     % instances are considered in the
                                     % evaluation of the sub. Otherwise,
                                     % unique (in terms of node sets) instances 
                                     % are taken into account (DEFAULT).
                                     % Also, redundancy is removed from the
                                     % main graph.
     options.subdue.supervised = false; % If true, graph search is performed over
				          % the whole data. If not, individual categories 
			                  % are searched, and the vocabularies are then 
			                  % combined.
end

