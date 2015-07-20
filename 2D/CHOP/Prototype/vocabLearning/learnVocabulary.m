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
%> @retval mainGraph The hierarchic object graphs.
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
function [ vocabulary, mainGraph, optimalThresholds, distanceMatrices, graphLevelIndices, edgeChangeLevel] = learnVocabulary( vocabLevel, graphLevel, leafNodes, ...
                                                            options, fileList)
    display('Vocabulary learning has started.');                          
    %% ========== Step 0: Set initial data structures ==========
    vocabulary = cell(options.maxLevels,1);
    mainGraph = cell(options.maxLevels,1);
    distanceMatrices = cell(options.maxLevels,1);
    graphLevelIndices = cell(options.maxLevels,1);
    optimalThresholds = single(repmat(options.subdue.threshold, options.maxLevels,1));
    
    %% ========== Step 1: Create first vocabulary and graph layers with existing node/edge info ==========
    %% Step 1.1: Prepare intermediate data structures for sequential processing.
    vocabulary(1) = {vocabLevel};
    mainGraph(1) = {graphLevel};
    
    %% Create distance matrices of the first level.
    if options.nodeSimilarityAllowed
        distanceMatrices(1) = {createDistanceMatrix(options.filters, options.distType, options.auto.deadFeatures)};
    else
        % Set dummy similarity matrix with only diagonals zero.
        numberOfNodes = numel(options.filters) - numel(options.auto.deadFeatures);
        newDistanceMatrix = inf(numberOfNodes, numberOfNodes, 'single');
        newDistanceMatrix(1:(numberOfNodes+1):numberOfNodes*numberOfNodes) = 0;
        distanceMatrices(1) = {newDistanceMatrix};
    end
    newDistanceMatrix = distanceMatrices{1};
    
    %% Get number of valid images in which we get gabor responses.
    numberOfImages = numel(unique([graphLevel.imageId]));
    
    % Open/close matlabpool to save memory. 
    if options.parallelProcessing
        s = matlabpool('size');
        if s>0
           matlabpool close; 
        end
        matlabpool('open', options.numberOfThreads);
    end
    
    %% Print first vocabulary and graph level.
    visualizeLevel( vocabulary{1}, [], [], [], [], 1, 0, 0, options);
    if ~isempty(distanceMatrices{1})
       imwrite(distanceMatrices{1}, [options.currentFolder '/debug/' options.datasetName '/level' num2str(1) '_dist.png']);
    end
    
    if options.debug
        matlabpool close; 
        display('........ Visualizing the realizations in the first level...');
        visualizeImages( fileList, vocabLevel, graphLevel, leafNodes, 1, options, 'train' );
        visualizeCroppedImgs( vocabulary{1}, 1, options);
        matlabpool('open', options.numberOfThreads);
    end
    printCloseFilters(distanceMatrices{1}, 1, options);
    
    %% Calculate statistics from this graph.
    display('........ Estimating statistics for level 1..');
    [avgShareability, avgCoverage] = saveStats(vocabLevel, graphLevel, leafNodes, numberOfImages, options, 'preInhibition', 1);
    display(['........ Average coverage of leaf nodes: ' num2str(avgCoverage) ', while average shareability is: ' num2str(avgShareability) ' percent.']); 
    
    %% Load categories. Analyzing categorization properties of the nodes. 
    load([options.currentFolder '/output/' options.datasetName '/export.mat']);
    firstLevelActivations = cat(1, mainGraph{1}.activation);
    
    % Read supervision flag for optimal threshold search.
    supervisedSelectionFlag = options.supervisedSelectionFlag;
    supervisedSelectionMode = options.supervisedSelectionMode;
    isSupervisedSelectionRunning = false;
    if supervisedSelectionFlag && strcmp(supervisedSelectionMode, 'manual')
        isSupervisedSelectionRunning = true;
    end
    previousAccuracy = 0;
    orgOptimizationFlag = options.optimizationFlag;
    
    %% ========== Step 2: Infer new parts by discovering frequent subs in data. ==========
    edgeChangeLevel = -1;
    for levelItr = 2:options.maxLevels
        % Obtain the pre-set threshold for this level, if there is one.
        if levelItr > (numel(options.subdue.presetThresholds) + 1)
            presetThreshold = options.subdue.threshold;
            options.optimizationFlag = orgOptimizationFlag;
        else
            presetThreshold = options.subdue.presetThresholds(levelItr-1);
            options.optimizationFlag = false;
        end
        
        %% Step 2.1: Run knowledge discovery to learn frequent compositions.
        [vocabLevel, graphLevel, optimalThreshold, isSupervisedSelectionRunning, previousAccuracy] = discoverSubs(vocabLevel, graphLevel, newDistanceMatrix, ...
            options, false, presetThreshold, levelItr-1, supervisedSelectionFlag, isSupervisedSelectionRunning, previousAccuracy);
        optimalThresholds(levelItr) = optimalThreshold;
        
        % Open/close matlabpool to save memory.
        matlabpool close;
        matlabpool('open', options.numberOfThreads);
        
        %% If no new subs have been found, finish processing.
        if isempty(vocabLevel)
           % Write previous level's appearances to the output folder.
           vocabulary = vocabulary(1:(levelItr-1),:);
           mainGraph = mainGraph(1:(levelItr-1),:);
           distanceMatrices = distanceMatrices(1:(levelItr-1),:);
           graphLevelIndices = graphLevelIndices(1:(levelItr-1),:);
           optimalThresholds = optimalThresholds(1:(levelItr-1),:);
           break; 
        end
        
        %% Assign realizations R of next graph level (l+1), and fill in their bookkeeping info.
        previousLevel = mainGraph{levelItr-1};
        graphLevel = fillBasicInfo(previousLevel, graphLevel, leafNodes, options.numberOfThreads);
        
        %% Calculate statistics from this graph.
        display('........ Before we apply inhibition, estimating statistics..');
        [avgShareability, avgCoverage] = saveStats(vocabLevel, graphLevel, leafNodes, numberOfImages, options, 'preInhibition', levelItr);
        display(['........ Average Coverage: ' num2str(avgCoverage) ', average shareability of compositions: ' num2str(avgShareability) ' percent.']); 
        
         %% If the subs have all been eliminated, finish processing.
        if isempty(vocabLevel)
           % Write previous level's appearances to the output folder.
           vocabulary = vocabulary(1:(levelItr-1),:);
           optimalThresholds = optimalThresholds(1:(levelItr-1),:);
           mainGraph = mainGraph(1:(levelItr-1),:);
           distanceMatrices = distanceMatrices(1:(levelItr-1),:);
           graphLevelIndices = graphLevelIndices(1:(levelItr-1),:);
           break; 
        end
        
        %% Inhibition! We process the current level to eliminate some of the nodes in the final graph.
        % The rules here are explained in the paper. Basically, each node
        % should introduce a novelty (cover a new area of the image). If it
        % fails to do so, the one which has a lower mdl ranking is
        % discarded. *Natural selection*
        if options.noveltyThr > 0
            display('........ Applying inhibition..');
            graphLevel=applyTestInhibition(graphLevel, options, levelItr);
        else
            display('........ No inhibition applied..');
        end
        % Open/close matlabpool to save memory.
        matlabpool close;
        matlabpool('open', options.numberOfThreads);
        
        %% Experimenting. After some point, we need to convert to centroid-based edge creation, no matter what.
        if avgCoverage < options.minContinuityCoverage && edgeChangeLevel == -1 && ~strcmp(options.edgeType, 'centroid')
            options.edgeType = 'centroid';
            display('........ Switching to -centroid- type edges!');
            edgeChangeLevel = levelItr;
        end
        
        %% Post-process graphLevel, vocabularyLevel to remove non-existent parts from vocabLevel.
        % In addition, we re-assign the node ids in graphLevel.
        if ~isempty(vocabLevel)
            display('........ Calculating distance matrix among the vocabulary nodes (in parallel)..');
            [vocabLevel, graphLevel, newDistanceMatrix, graphLabelAssgnArr] = postProcessParts(vocabLevel, graphLevel, distanceMatrices{levelItr-1}, options);
            graphLevelIndices{levelItr} = graphLabelAssgnArr;
            distanceMatrices{levelItr} = newDistanceMatrix;
        else
            newDistanceMatrix = [];
        end
        
        %% Calculate statistics from this graph.
        display('........ Estimating post-inhibition statistics..');
        [avgShareability, avgCoverage] = saveStats(vocabLevel, graphLevel, leafNodes, numberOfImages, options, 'postInhibition', levelItr);
        
        % display debugging info.
        display(['........ Inhibition applied with novelty thr: ' num2str(options.noveltyThr) ' and edge novelty thr: ' num2str(options.edgeNoveltyThr) '.']);
        display(['........ Remaining: ' num2str(numel(graphLevel)) ' realizations belonging to ' num2str(numel(vocabLevel)) ' compositions.']);
        display(['........ Average Coverage: ' num2str(avgCoverage) ', average shareability of compositions: ' num2str(avgShareability) ' percent.']); 
        
        %% Step 2.4: In order to do proper visualization, we learn precise positionings of children for every vocabulary node.
        vocabLevel = learnChildPositions(vocabLevel, graphLevel, previousLevel);
        
        %% Step 2.5: Create the parent relationships between current level and previous level.
        vocabulary = mergeIntoGraph(vocabulary, vocabLevel, leafNodes, levelItr, 0);
        mainGraph = mergeIntoGraph(mainGraph, graphLevel, leafNodes, levelItr, 1);
        
        if levelItr == options.maxLevels
            vocabulary = vocabulary(1:(levelItr),:);
            optimalThresholds = optimalThresholds(1:(levelItr),:);
            mainGraph = mainGraph(1:(levelItr),:);
            distanceMatrices = distanceMatrices(1:(levelItr),:);
            graphLevelIndices = graphLevelIndices(1:(levelItr),:);
            break;
        end
        
        %% Step 2.6: Create object graphs G_(l+1) for the next level, l+1.
        % Extract the edges between new realizations to form the new object graphs.
        [mainGraph] = extractEdges(mainGraph, options, levelItr);
        graphLevel = mainGraph{levelItr};
        
        %% Print vocabulary and graph level to output images (reconstruction).
        if ~isempty(newDistanceMatrix)
           imwrite(newDistanceMatrix, [options.currentFolder '/debug/' options.datasetName '/level' num2str(levelItr) '_dist.png']);
        end
        if ~isempty(vocabLevel)
            display('........ Visualizing previous levels...');
            visualizeLevel( vocabLevel, vocabulary, graphLevel, firstLevelActivations, leafNodes, levelItr, numel(vocabulary{1}), numel(vocabulary{levelItr-1}), options);
            if options.debug
               display('........ Visualizing realizations on images...');
               if ~isempty(vocabLevel)
                    matlabpool close;
                    visualizeImages( fileList, vocabLevel, graphLevel, leafNodes, levelItr, options, 'train' );
                    visualizeCroppedImgs( vocabLevel, levelItr, options);
                    matlabpool('open', options.numberOfThreads);
               end
            end
            printCloseFilters(newDistanceMatrix, levelItr, options); 
        end
        
        % Open/close matlabpool to save memory.
        matlabpool close;
        matlabpool('open', options.numberOfThreads);
        
        %% Step 2.6: If no new edges found, kill program.
        newEdgesAvailable = ~isempty(cat(1, mainGraph{levelItr}.adjInfo));
        if ~newEdgesAvailable
            vocabulary = vocabulary(1:(levelItr),:);
            optimalThresholds = optimalThresholds(1:(levelItr),:);
            mainGraph = mainGraph(1:(levelItr),:);
            distanceMatrices = distanceMatrices(1:(levelItr),:);
            graphLevelIndices = graphLevelIndices(1:(levelItr),:);
        end
    end
end