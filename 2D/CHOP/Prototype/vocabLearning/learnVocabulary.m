%> Name: learnVocabulary
%>
%> Description: Given the graph of first-level responses and their
%> relations, this function extracts a hierarchical vocabulary by detecting
%> structural patterns inherent in the input graph. The pattern discovery
%> phase is hierarchic, and continues until a the graph is compressed into
%> a single node.
%>
%> @param vocabLevel The first level of vocabulary, level 1 parts.
%> @param graphLevel Object graphs level 1.
%> @param options Program options
%> @param fileList input image name list.
%>
%> @retval vocabulary The hierarchic vocabulary learnt from the data.
%> @retval distanceMatrices Distance matrix of compositions.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 26.11.2013
%> Ver 1.1 on 02.12.2013 Completed unlimited number of graph generation.
%> Ver 1.2 on 12.01.2014 Commentary changes, for unified code look.
%> Ver 1.3 on 28.01.2014 Mode calculation put in a separate file.
%> Ver 1.4 on 03.02.2014 Refactoring
function [] = learnVocabulary( vocabLevel, vocabLevelDistributions, graphLevel, leafNodes, leafNodeCoords, allResponses, ...
                                                            options, fileList)
    display('Vocabulary learning has started.');                          
    %% ========== Step 0: Set initial data structures ==========
    vocabulary = cell(options.maxLevels,1);
    vocabularyDistributions = cell(options.maxLevels,1);
    allModes = cell(options.maxLevels,1);
    distanceMatrices = cell(options.maxLevels,1);
    modeProbs = cell(options.maxLevels,1);
    
    %% ========== Step 1: Create first vocabulary and graph layers with existing node/edge info ==========
    %% Step 1.1: Prepare intermediate data structures for sequential processing.
    vocabulary(1) = {vocabLevel};
    vocabularyDistributions(1) = {vocabLevelDistributions};
 %   allModes(1) = {modes};
 %   modeProbs(1) = {modeProbArr};
    
    %% Create distance matrices of the first level.
    [eucDistanceMatrix, nodeDistanceMatrix] = createDistanceMatrix(options.filters, options.filterType);
    distanceMatrices(1) = {nodeDistanceMatrix};
    newDistanceMatrix = distanceMatrices{1};
    
    %% Get number of valid images in which we get gabor responses.
    numberOfImages = numel(unique([graphLevel.imageId]));
    
    % Open threads for parallel processing.
    if options.parallelProcessing && usejava('jvm')
        if exist('parpool', 'file')
            p = gcp('nocreate');
            if ~isempty(p)
                delete(p);
            end
            parpool(options.numberOfThreads);
        else
            s = matlabpool('size');
            if s>0
               matlabpool close; 
            end
            matlabpool('open', options.numberOfThreads);
        end
    end
    
    %% Print first vocabulary and graph level.
    allNodeInstances = visualizeLevel( vocabulary{1}, [], [], [], [], 1, 0, 0, options);
    if ~isempty(distanceMatrices{1})
       imwrite(distanceMatrices{1}, [options.currentFolder '/debug/' options.datasetName '/level' num2str(1) '_dist.png']);
    end
    
    if options.debug
%        matlabpool close; 
        display('........ Visualizing the realizations in the first level...');
        visualizeImages( fileList, graphLevel, [], allNodeInstances, leafNodes, leafNodeCoords, 1, options, 'train' );
%        visualizeCroppedImgs( vocabulary{1}, 1, options);
%        matlabpool('open', options.numberOfThreads);
    end
    printCloseFilters(eucDistanceMatrix, [], 1, options);
    
    %% Calculate statistics from this graph.
    display('........ Estimating statistics for level 1..');
    [avgShareability, avgCoverage, maxCoverageVals] = saveStats(vocabLevel, graphLevel, leafNodeCoords, [], numberOfImages, options, 'preInhibition', 1);
    display(['........ Average coverage of leaf nodes: ' num2str(avgCoverage) ', while average shareability is: ' num2str(avgShareability) ' percent.']);
    
    %% Load categories. Analyzing categorization properties of the nodes. 
    load([options.currentFolder '/output/' options.datasetName '/export.mat']);
    firstLevelActivations = cat(1, graphLevel.activation);
    
    % Read supervision flag for optimal threshold search.
    supervisedSelectionFlag = options.supervisedSelectionFlag;
    supervisedSelectionMode = options.supervisedSelectionMode;
    isSupervisedSelectionRunning = false;
    if supervisedSelectionFlag && strcmp(supervisedSelectionMode, 'manual')
        isSupervisedSelectionRunning = true;
    end
    previousAccuracy = 0; %#ok<*NASGU>
    
    %% Obtain coords.
    level1Coords = [cat(1, graphLevel.imageId), int32(cat(1, graphLevel.precisePosition)), cat(1, graphLevel.position)];
    firstLevelPrecisePositions = cat(1, graphLevel.precisePosition);
    
    %% If we have passed the limit for number of images, we go into a different mode 
    % where we perform matlab restarts for memory optimization.
    if numel(fileList) >= options.matlabRestartImageCount
         options.restartFlag = true;
    else
         options.restartFlag = false;
    end
    
    %% Set number of category level OR nodes.
    options.reconstruction.numberOfORNodesCategoryLayer = options.articulationsPerCategory * options.numberOfCategories;
    
    %% Export first level realizations.
    levelItr = 1;
    disp('Exporting first level realizations.');
    [preExportArr, preActivationArr, prePooledPositions, preLeafNodes] = exportRealizations(graphLevel, levelItr); %#ok<ASGLU>
      
    %% ========== Step 2: Infer new parts by discovering frequent subs in data. ==========
    edgeChangeLevel = -1;
    bInfo = dbstatus;
    
    % Set the initial stop flag to 0. Algorithm will automatically stop
    % when most objects are covered.
    options.stopFlag = false;
    datasetName = options.datasetName;
    
    % Learn stride.
    if strcmp(options.filterType, 'gabor')
         stride = options.gabor.stride;
    else
         stride = options.auto.stride;
    end
    poolDim = options.poolDim;
    
    % Save workspace into a file.
    if numel(fileList) >= options.saveWorkspaceImageCount
         disp('Saving layer 1 workspace.');
         saveWorkspace();
    end
    
    if options.restartFlag
         
         % Reduce memory consumption by writing all stuff to files,
         % clearing all, and then returning back to computation.
         scriptName = [options.currentFolder '/' options.datasetName '.sh'];
         
         % Let's generate the restart script.
         fid = fopen(scriptName, 'w');
         if ismac
              versionName = version('-release');
              matlabString = ['export PATH=/Applications/MATLAB_R' versionName '.app/bin:$PATH'];
         else
              matlabString = [];
         end
         
         % If no jvm is requested, move on.
         if ~usejava('jvm')
             argString = '-nosplash -nojvm -nodisplay';
         else
             argString = [];
         end
         
         fileString = fprintf(fid, '#!/bin/bash\n%s\nscreen -d -m matlab -r %s "learnVocabularyLevel()"', matlabString, argString);
         system(['chmod 755 ' scriptName]);
         fclose(fid);
         system(scriptName);
         exit
    else
         learnVocabularyLevel();
    end
end