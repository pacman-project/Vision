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
    if ~isTraining
        options.debug = true;
    else
        options.debug = true;           % If debug = 1, additional output will be 
                                 % generated to aid debugging process.
    end
    options.fastStatLearning = false;
    options.vis.printTrainRealizations = true;
    options.backgroundClass = 'Background'; % The string that identifies 
                                            % background class. Images from
                                            % this set will be used as
                                            % negative examples in
                                            % training.
    %% ========== TO BE MIGRATED TO INDIVIDUAL PARAMETER FILES ==========                    
    options.autoNormalize = 'whiten'; %Or 'normalize'
    options.poolDim = 2; % The pooling window size. Could be 1, 2, 3...
                                      % If 2, size of the receptive field
                                      % is halved, as in the previous
                                      % version of the algorithm.
 %   options.noPoolingLayers = [4 6 8 10 12 14]; % No pooling is performed over the outcome of these layers.
%    options.smallRFLayers = [3 5 7 9 11 13];
    options.noPoolingLayers = [3 5 7 9 10 11 12];
    options.smallRFLayers = [2 4 6 8];
    options.poolFlag = true; % If true, we apply max pooling to reduce number of 
                                            %  realizations. Otherwise,
                                            %  positions are updated, but
                                            %  no elimination is performed.
    
    %% ========== PART SELECTION PARAMETERS ==========
    options.partSelectionFlag = true; % When this flag is false, the default 
                                      % part selection scheme is applied
                                      % (MDL-based). If true, a second pass
                                      % (either based on reconstruction or
                                      % supervision) is applied to the
                                      % shortlisted part candidates that
                                      % are obtained using SUBDUE's default
                                      % part selection mechanism.
    options.optimizationFlag = false; % If this flag is true, threshold 
                                     % optimization will take place at the 
                                     % end of SUBDUE-based part search. 
    options.supervisedSelectionFlag = false; % If true, the algorithm will 
                                   % switch to supervised part selection.
    options.supervisedSelectionMode = 'manual'; % If 'auto', the system will 
                                   % switch to discriminative threshold search 
                                   % after first performance drop in
                                   % unsupervised learning. If 'manual', it
                                   % will go supervised since level 2.
    options.subdue.presetThresholds = []; % To be set later in parameter files.
    
    %% ========== GRAPH MATCHING PARAMETERS ==========
    options.nodeMatchingFlag = 0; % If 0, only nodes with the same node ids are matched. 
                                  % If 1, elastic node matching is
                                  % performed.
    options.edgeMatchingFlag = 1; % If 0, only edges with the same id can be
                                  % matched. If 1, elastic edge matching is
                                  % performed. 
                                   
    %% ========== PARALLEL PROCESSING PARAMETERS ==========
    options.parallelProcessing = true;
    options.numberOfThreads = min(feature('NumCores')*2-1,12);
%    options.numberOfThreads = 12;
    if exist('parpool', 'file')
        options.parpoolFlag = true;
    else
        options.parpoolFlag = false;
    end
        
    
    %% ========== DATASET-SPECIFIC PROGRAM PARAMETERS ==========
    % Learn dataset path relative to this m file
    currentFileName = mfilename('fullpath');
    [currentPath, ~, ~] = fileparts(currentFileName);
    cd([currentPath '/parameters/']);
    if exist(['SetParameters' datasetName '.m'], 'file')
        fh = str2func(['SetParameters' datasetName]);
    else
        fh = @SetParametersCommon;
    end
    options = fh(datasetName, options);
    cd(currentPath);
    
    %% ========== FOLDER STRUCTURE INITIALIZATION ==========
    % Set folder parameters.
    options.currentFolder = currentPath;
    options.debugFolder = [currentPath '/debug/' datasetName];
    options.processedFolder = [currentPath '/output/' datasetName '/original'];
    options.processedGTFolder = [currentPath '/output/' datasetName '/gt'];
    options.outputFolder = [currentPath '/output/' datasetName];
    options.testOutputFolder = [options.outputFolder '/test'];
    options.testInferenceFolder = [options.outputFolder '/test/inference'];
    options.smoothedFolder = [currentPath '/output/' datasetName '/smoothed'];
    options.CNNFolder = [currentPath '/CNN/' datasetName];
    
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
    options.graphNode = GraphNode;
    
    %% ========== LOW - LEVEL FILTER GENERATION ==========
    [filters, filterImages] = createFilters(options);
    options.filters = filters;
    options.filterImages = filterImages;
    options.numberOfFilters = numel(filters);
    
    %%  ========== SINGLE PRECISION  ========== 
    options.singlePrecision = single(0.0001);
    
    %% ========== FILTER MATRIX & DATA STRUCTURES GENERATION ==========
    % We create a feature matrix out of these filters for fast
    % processing.
    filterMatrix = zeros(options.numberOfFilters, numel(filters{1}));
    stDevs = zeros(options.numberOfFilters,1);
    if size(filterMatrix,2) >0
        for filtItr = 1:options.numberOfFilters
            filter1 = filters{filtItr};
            filterMatrix(filtItr,:) = filter1(:);
            stDev = 0;
            for bandItr = 1:size(filter1,3)
                bandImg = filter1(:,:,bandItr);
                stDev = stDev + std(bandImg(:));
            end
            stDevs(filtItr) = stDev / size(filter1,3);
        end
    end
    options.filterMatrix = filterMatrix;
    stDevs = stDevs / max(stDevs);
    
    if strcmp(options.filterType, 'auto') 
        if exist([options.currentFolder '/filters/vis/' datasetName '/C.mat'], 'file') 
            load([options.currentFolder '/filters/vis/' datasetName '/C.mat'], 'whMat', 'mu', 'invMat');
            options.auto.whMat = whMat;
            options.auto.invMat = invMat;
            options.auto.mu = mu;

            % Mark dead features.
            options.auto.deadFeatures = find(stDevs < options.auto.deadFeatureStd );
        else
            display('Auto-features not learned. Please run automatic feature learning before running CHOP. You can disregard this message if you are running that already.');
        end
    else
        options.auto.deadFeatures = [];
    end
end

