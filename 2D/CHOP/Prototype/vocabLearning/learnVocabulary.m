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
%> @param modes The mode array that only includes first level relations.
%> @param options Program options
%> @param fileList input image name list.
%>
%> @retval vocabulary The hierarchic vocabulary learnt from the data.
%> @retval redundantVocabulary The vocabulary which contains redundant
%> compositions, which have indices of actual compositions in their 'label'
%> fields.
%> @retval mainGraph The hierarchic object graphs.
%> @retval modes Cell array including modes belonging to each layer, each 
%>      in a separate cell.
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
function [ vocabulary, redundantVocabulary, mainGraph, modes, distanceMatrices] = learnVocabulary( vocabLevel, graphLevel, leafNodes, modes, ...
                                                            options, fileList)
    display('Vocabulary learning has started.');                          
    %% ========== Step 0: Set initial data structures ==========
    vocabulary = cell(options.maxLevels,1);
    redundantVocabulary = cell(options.maxLevels,1);
    mainGraph = cell(options.maxLevels,1);
    distanceMatrices = cell(options.maxLevels,1);
    
    %% ========== Step 1: Create first vocabulary and graph layers with existing node/edge info ==========
    %% Step 1.1: Prepare intermediate data structures for sequential processing.
    vocabulary(1) = {vocabLevel};
    mainGraph(1) = {graphLevel};
    previousModes = [];
    
    %% Create distance matrices of the first level.
    distanceMatrices(1) = {createDistanceMatrix(options.filters)};
    
    %% Get number of valid images in which we get gabor responses.
    numberOfImages = numel(unique([graphLevel.imageId]));
    
    %% Print first vocabulary and graph level.
    visualizeLevel( vocabulary{1}, [], [], [], 1, previousModes, 0, options, 0);
    if options.debug
        display('........ Visualizing the realizations in the first level...');
        visualizeImages( fileList, vocabLevel, graphLevel, leafNodes, 1, options, 'train' );
    end
    
    %% Calculate statistics from this graph.
    display('........ Estimating statistics for level 1..');
    [avgShareability, avgCoverage] = saveStats(vocabLevel, graphLevel, leafNodes, numberOfImages, options, 'preInhibition', 1);
    display(['........ Average coverage of leaf nodes: ' num2str(avgCoverage) ', while average shareability is: ' num2str(avgShareability) ' percent.']); 
    
    %% ========== Step 2: Infer new parts by discovering frequent subs in data. ==========
    for levelItr = 2:options.maxLevels
        %% Step 2.0: Assign current modes.
        currentModes = [];
        if ~isempty(modes) && numel(modes)>=(levelItr-1)
            currentModes = modes{levelItr-1};
        end
        
        if strcmp(options.subdue.evalMetric, 'mdl')
            allEdges = {graphLevel.adjInfo}';
            graphSize = sum(cellfun(@(x) size(x,1), allEdges)) * options.subdue.mdlEdgeWeight + ...
                numel(graphLevel) * options.subdue.mdlNodeWeight;
            clear allEdges;
        else
            graphSize = 1;
        end
        
        %% Step 2.1: Run knowledge discovery to learn frequent compositions.
        [vocabLevel, graphLevel] = discoverSubs(vocabLevel, redundantVocabulary{levelItr-1}, graphLevel, ...
            options, false, levelItr-1);
        
        %% If no new subs have been found, finish processing.
        if isempty(vocabLevel)
           % Write previous level's appearances to the output folder.
           vocabulary = vocabulary(1:(levelItr-1),:);
           redundantVocabulary = redundantVocabulary(1:(levelItr-1),:);
           mainGraph = mainGraph(1:(levelItr-1),:);
           distanceMatrices = distanceMatrices(1:(levelItr-1),:);
           break; 
        end
        
        %% Assign realizations R of next graph level (l+1), and fill in their bookkeeping info.
        previousLevel = mainGraph{levelItr-1};
        graphLevel = fillBasicInfo(previousLevel, graphLevel, leafNodes);
        
        %% Calculate statistics from this graph.
        display('........ Before we apply inhibition, estimating statistics..');
        [avgShareability, avgCoverage] = saveStats(vocabLevel, graphLevel, leafNodes, numberOfImages, options, 'preInhibition', levelItr);
        display(['........ Average Coverage: ' num2str(avgCoverage) ', average shareability of compositions: ' num2str(avgShareability) ' percent.']); 
        
        %% Combining parts to have more generic parts. 
        % Apologies for this section to be 'hidden' for now. It'll be unhid as soon as possible.
        display('........ Combining parts..');
        [vocabLevel, graphLevel, newDistanceMatrix, subClasses] = combineParts(vocabLevel, graphLevel, currentModes, distanceMatrices{levelItr-1}, graphSize, options);
        
        %% Inhibition! We process the current level to eliminate some of the nodes in the final graph.
        % The rules here are explained in the paper. Basically, each node
        % should introduce a novelty (cover a new area of the image). If it
        % fails to do so, the one which has a lower mdl ranking is
        % discarded. *Natural selection*
        display('........ Applying inhibition..');
        graphLevel=applyTestInhibition(graphLevel, options, levelItr);
        
        %% Post-process graphLevel, vocabularyLevel and find redundantVocabularyLevel.
        % This part is again related to combining parts.
        % Apologies for this section to be 'hidden' for now. It'll be unhid as soon as possible.
        [vocabLevel, redundantVocabLevel, graphLevel, newDistanceMatrix] = postProcessParts(vocabLevel, subClasses, newDistanceMatrix, graphLevel);
        distanceMatrices{levelItr} = newDistanceMatrix;
        
        %% Calculate statistics from this graph.
        [avgShareability, avgCoverage] = saveStats(vocabLevel, graphLevel, leafNodes, numberOfImages, options, 'postInhibition', levelItr);
        
        % display debugging info.
        display(['........ Inhibition applied with novelty thr: ' num2str(options.noveltyThr) ' and edge novelty thr: ' num2str(options.edgeNoveltyThr) '.']);
        display(['........ Remaining: ' num2str(numel(graphLevel)) ' realizations belonging to ' num2str(numel(vocabLevel)) ' compositions.']);
        display(['........ Average Coverage: ' num2str(avgCoverage) ', average shareability of compositions: ' num2str(avgShareability) ' percent.']); 
        
        %% Step 2.4: Create the parent relationships between current level and previous level.
        vocabulary = mergeIntoGraph(vocabulary, vocabLevel, leafNodes, levelItr, 0);
        redundantVocabulary(levelItr) = {redundantVocabLevel};
        mainGraph = mergeIntoGraph(mainGraph, graphLevel, leafNodes, levelItr, 1);
        
        %% Step 2.5: Create object graphs G_(l+1) for the next level, l+1.
        % Extract the edges between new realizations to form the new object graphs.
        [modes, mainGraph] = extractEdges(mainGraph, options, levelItr, modes);
        graphLevel = mainGraph{levelItr};
        
        %% Print vocabulary and graph level to output images (reconstruction).
        if ~isempty(modes)
            display('........ Visualizing previous level...');
            visualizeLevel( vocabLevel, graphLevel, leafNodes, distanceMatrices{1}, levelItr, modes{levelItr-1}, numel(vocabulary{1}), options, 0);
        end
        if options.debug
           display('........ Visualizing realizations on images...');
           if ~isempty(vocabLevel)
               visualizeImages( fileList, vocabLevel, graphLevel, leafNodes, levelItr, options, 'train' );
           end
        end
        
        %% Step 2.6: If no new edges found, kill program.
        newEdgesAvailable = ~isempty(cat(1, mainGraph{levelItr}.adjInfo));
        if ~newEdgesAvailable
            vocabulary = vocabulary(1:(levelItr),:);
            redundantVocabulary = redundantVocabulary(1:(levelItr),:);
            mainGraph = mainGraph(1:(levelItr),:);
            distanceMatrices = distanceMatrices(1:(levelItr),:);
        end
    end
end