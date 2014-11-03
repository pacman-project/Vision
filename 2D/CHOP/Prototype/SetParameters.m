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
        options.debug = 1;
    else
        options.debug = 1;           % If debug = 1, additional output will be 
                                 % generated to aid debugging process.
    end                             
    options.backgroundClass = 'Background'; % The string that identifies 
                                            % background class. Images from
                                            % this set will be used as
                                            % negative examples in
                                            % training.
                                 
    %% ========== PARALLEL PROCESSING PARAMETERS ==========
    options.parallelProcessing = true;
    options.numberOfThreads = min(feature('NumCores'),12); % For my macbook pro
    
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
    
    %% ========== FILTER MATRIX & DATA STRUCTURES GENERATION ==========
    if strcmp(options.filterType, 'auto') 
        if exist([options.currentFolder '/filters/vis/' datasetName '/C.mat'], 'file') 
            load([options.currentFolder '/filters/vis/' datasetName '/C.mat'], 'whMat', 'mu', 'invMat');
            options.auto.whMat = whMat;
            options.auto.invMat = invMat;
            options.auto.mu = mu;

            % We create a feature matrix out of these filters for fast
            % processing.
            filterMatrix = zeros(options.numberOfFilters, numel(filters{1}));
            if size(filterMatrix,2) >0
                for filtItr = 1:options.numberOfFilters
                    filterMatrix(filtItr,:) = filters{filtItr}(:);
                end
            end
            options.filterMatrix = filterMatrix;

            stDevs = std(filterMatrix, 0, 2);
            stDevs = stDevs / max(stDevs);
            % Mark dead features.
            options.auto.deadFeatures = find(stDevs < options.auto.deadFeatureStd );
        else
            display('Auto-features not learned. Please run automatic feature learning before running CHOP. You can disregard this message if you are running that already.');
        end
    end
end

