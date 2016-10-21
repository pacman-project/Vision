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
function [ options ] = SetParametersGWP_2k( datasetName, options )

    %% ========== DATASET - RELATED PARAMETERS ==========
    options.datasetName = datasetName;
    options.learnVocabulary = 1; % If 1, new vocabulary is learned. 
    options.testImages = 1;      % If 1, the test images are processed.
    options.numberOfGaborFilters = 8; % Number of Gabor filters at level 1.
    options.matlabRestartImageCount = 100000; % If we have more than this number 
    % of images, matlab restarting is performed after every layer.
    
        %% ========== LOW - LEVEL FILTER PARAMETERS ==========
    options.filterType = 'gabor'; % If 'gabor': Steerable Gabor filters used 
                                  % as feature detectors.
                                  % If 'auto': Autodetected features.
                                  % Random patches are clustered to obtain
                                  % a number of unsupervised features.
    options.gaborFilterThr = 0.1; % Min response threshold for convolved features, 
                                  % taken as the percentage of max response 
                                  % in each image.
    options.innerGaborFilterThr = 0.05; % Min response threshold for inner contours.
    options.absGaborFilterThr = 0; % Absolute response threshold for low-level  % For depth!
                                   % responses. ~80 for natural images 
                                   % (depends on many factors though, including 
                                   % size of the filter).
    options.gaborFilterSize = 13;       % Size of a gabor filter. Please note 
                                        % that the size also depends on the 
                                        % filter parameters, so consider them 
                                        % all when you change this!
%    options.gabor.sigma = 2.5;            % Gabor filter parameters
 %   options.gabor.sigma = 1.21;    
    options.gabor.sigma = 1;    
    options.gabor.theta = 0;
    options.gabor.lambda = 1;
    options.gabor.psi = 0;
    options.gabor.gamma = 0.2;
    options.gabor.inhibitionRadius = 2;
                                        % The inhibition radius basically 
                                        % defines the half of the square's
                                        % size in which weaker responses other 
                                        % than the seed node will
                                        % be surpressed. The image is
                                        % pooled with a stride of
                                        % (inhibitionRadius+1) after
                                        % inhibition.
    options.gabor.stride = 2;           % Stride to use when extracting gabor
                                       % features.     
    options.autoFilterSize = 10;         % Size (one side) of a autodetected 
                                        % filter. Assumed to be NxNxD.
    options.auto.inhibitionRadius = 2;
    options.autoFilterThr = 0;       % Min response threshold for convolved 
                                       % features, assigned as this percentage 
                                       % of the max response in each image.
    options.autoFilterCount = 100;      % Number of auto-detected filters.
    options.autoFilterPatchCount = 10000; % Number of random patches used 
                                           % to find auto-detected filters.
    options.auto.stride = 3;           % Stride to use when extracting first-
                                       % level features. Only works in
                                       % auto-filter mode, since gabors are
                                       % extracted using conv2, convolution
                                       % implementation of matlab.d                                 
    options.auto.deadFeatureStd = 0; % In case of auto-learned features, 
                                       % some dead features may come up.
                                       % The standard deviation check is
                                       % used to eliminate uniform
                                       % features, assigned as this percentage 
                                       % of the max std dev in filters.
    options.distType = 'sift'; % Either 'modal', 'hog' or 'hu', or 'sift'.
    if strcmp(options.filterType, 'gabor')
         options.stride = options.gabor.stride;
    else
         options.stride = options.auto.stride;
    end
    
    %% ========== GT Parameters ==========
    options.useGT = true;              % If true, gt info is used. 
    options.gtType = 'contour';        % 'contour' type gt: nodes lying under
                                       % the gt contour is examined (within
                                       % a neighborhood defined by
                                       % contourGTNeighborhood). 
                                       % 'bbox' type gt: nodes in the
                                       % gt bounding box are examined.
    options.contourGTNeighborhood = 0;% width of the band along the contour 
                                       % (half width, actual width is
                                       % double this value)  in which nodes
                                       % and relations are examined.

    %% ========== CRUCIAL METHOD PARAMETERS (COMPLEXITY, RELATIONS) ==========
    options.noveltyThr = 0;           % The novelty threshold used in the 
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
    options.maxShareability = 0.2;       % If more than this percent of a nodes' 
                                             % leaf nodes are shared, this node is discarded.
    options.minEdgeNoveltyThr = 0.0;    % Minimum edge novelty threshold.
    options.edgeNoveltyThrRate = 0.0;   % The edge novelty threshold 
                                       % is supposed to reduce each level by 
                                       % this amount.
    options.edgeType = 'continuity';     % If 'centroid', downsampling is
                                       % applied at each layer, and edges
                                       % link spatially adjacent (within
                                       % its neighborhood) nodes.
                                       % If 'continuity', linking depends
                                       % on another condition: Leaf node
                                       % sets of two nodes should be linked
                                       % by at least one edge in the first
                                       % level graph. This option can be
                                       % used in the first few layers in
                                       % order to ensure continuity (e.g.
                                       % smooth boundaries/surfaces) in upper 
                                       % layers. 
    options.minContinuityCoverage = 0.94; % If data coverage drops below this,
                                         % we switch to 'centroid' nodes.
    options.missingNodeThr = 0; % Each node should cover this percentage of the nodes in its RF.
    options.categoryLevelMissingNodeThr = 0; %In category level, missing node threshold is treated differently.
    options.edgeChangeLevel = -1; % Going to be updated in the code later.
    options.maxEdgeChangeLevel = 10; % If this is the layer we're working on, we switch to centroid edges.
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
        options.receptiveFieldSize = 5; % DEFAULT 5
    else
        options.receptiveFieldSize = 13;
    end                                  % Size (one side) of the receptive field at
                                         % first level. Please note that in
                                         % each level of the hierarchy, the
                                         % receptive field size grows by 
                                         % 1/scaling.
                                         
    options.smallReceptiveFieldSize = 9; % Small receptive field size that is used in some levels.
    options.circularRF = true; % If true, the RF is treated circularly. Otherwise it's square.
    options.maxNodeDegree = 20;        % (N) closest N nodes are linked for 
                                       % every node in the object graphs.
    options.maxFirstLevelNodeDegree = 10;  % (N) closest N nodes are linked for 
                                       % every node in the object graphs (layer 1).
    options.minimalEdgeCount = true; % If true, coverage-based edge creation 
                                                            % is performed. If an edge is not contributing to 
                                                            % coverage, it's not generated.
    options.maxImageDim = 2000; %Max dimension of the 
                                       % images the algorithm will work
                                       % with. If one size of a image in
                                       % the dataset is larger than this
                                       % value, it will be rescaled to fit
                                       % in a square of
                                       % maxImageDim x maxImageDim. Aspect ratio
                                       % will be preserved. Set to a large
                                       % value to avoid rescaling.
    options.maxLevels = 8;    % The maximum level count for training.
    options.maxInferenceLevels = 8; % The maximum level count for testing.
    
    %% ========== INFERENCE PARAMETERS ==========
    options.fastInference = true; % If set, faster inference (involves 
                                  % inhibition) is performed.
    options.favorParam = 1;      % Between 1:100, if it increases, category 
                                 % nodes with peaks towards a single category 
                                 % will be favoured more, rather than those 
                                 % with relatively uniform distribution.
                                 % Used in determining the category of a node.
                                 
    options.categoryLevel = 8; % The level where we switch from 
                                                  % geometry-based grouping
                                                  % to category nodes.
                                                  % In effect, the number
                                                  % of nodes in the
                                                  % vocabulary is reduced
                                                  % drastically. Set to an
                                                  % unlikely value (100)
                                                  % for no category nodes.
    options.categoryLevelCoverage = 0.8; % If each node is covering this percent 
                                                                  % of its associated image, we switch to 
                                                                  % category node grouping.
    options.articulationsPerCategory = 30; % Number of top level nodes 
               % for each category. Category nodes will be formed by
               % grouping top level parts.
               % Hint: A good number is the number of poses per object.

    %% ========== RECONSTRUCTION PARAMETERS ==========
    options.reconstruction.numberOfReconstructiveSubs = 5000; % The maximum 
                                           % number of reconstructive parts
                                           % that can be selected.
%    options.reconstruction.numberOfORNodes = [100 120 250 500 750 1000 1000 1000 1000 1000 1000 1000 1000 1000]; % The maximum 
                                           % number of reconstructive parts
                                           % that can be selected.
    options.reconstruction.numberOfORNodes = []; % The maximum 
                                           % number of reconstructive parts
                                           % that can be selected.
    options.reconstruction.maxNumberOfORNodes = 3000; % Ideal number of OR 
    % nodes is searched in the range of
    % minNumberOfOrNodes-maxNumberOfOrNodes. If cannot be found,
    % numberOfORNodes (set above) is used.
    
    options.reconstruction.minNumberOfORNodes = 50; % If the number of 
                                   % nodes is less than this value, we don't need compression.
                                           
    %% ========== GRAPH MATCHING PARAMETERS ==========
    options.nodeSimilarityAllowed = false; % If true, node similarities are 
                                           % considered in graph matching.
                                           % If not, identicality in labels
                                           % represents zero-cost matching,
                                           % while every other kind of node
                                           % correspondance yields a cost
                                           % of 1 (max value). 
    options.edgeSimilarityAllowed = false;  % If true, edge similarities are 
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
    options.subdue.mdlNodeWeight = 12;      % Weight of a node in DL calculations 
                                           % in MDL-based evaluation
                                           % metric. Cost of a node =
                                           % labelId (int, 4 byte) + pointer to
                                           % edges (int, 4 byte) = 8.
    options.subdue.mdlEdgeWeight = 8;      % Weight of an edge in DL calculations 
                                           % in MDL-based evaluation
                                           % metric. Cost of an edge =
                                           % edgeLabelId (int, 4 byte) + 
                                           % destinationNode (int,4 byte) + 
                                           % isDirected (byte, 1 byte) = 9.
    options.subdue.maxTime = 1800;          % Max. number of seconds subdue is
                                            % allowed to run. Typically
                                            % around 100 (secs) for toy data. 
                                            % You can set to higher values
                                            % (e.g. 3600 secs) for large
                                            % datasets.
    options.subdue.minSize = 1; % Minimum number of nodes in a composition.
    options.subdue.maxSize = min(4, (options.maxNodeDegree + 1)); % Maximum number of nodes in a composition.
    options.subdue.nsubs = 50000;  % Maximum number of nodes allowed in a level. Subs of size 1, 2, ..., n each can contribute in this number.
    options.subdue.beam = 5000;   % Beam length in SUBDUE' search mechanism.
    options.subdue.overlap = false;   % If true, overlaps between a substructure's 
                                     % instances are considered in the
                                     % evaluation of the sub. Otherwise,
                                     % unique (in terms of node sets) instances 
                                     % are taken into account (DEFAULT).
     options.subdue.supervised = false; % If true, graph search is performed over
				          % the whole data. If not, individual categories 
			                  % are searched, and the vocabularies are then 
			                  % combined.
end

