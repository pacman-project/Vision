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
%> @retval mainGraph The hierarchic object graphs.
%> @retval modes Cell array including modes belonging to each layer, each 
%>      in a separate cell.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 26.11.2013
%> Ver 1.1 on 02.12.2013 Completed unlimited number of graph generation.
%> Ver 1.2 on 12.01.2014 Commentary changes, for unified code look.
%> Ver 1.3 on 28.01.2014 Mode calculation put in a separate file.
%> Ver 1.4 on 03.02.2014 Refactoring
function [ vocabulary, mainGraph, modes, allOppositeModes, highLevelModes ] = learnVocabulary( vocabLevel, graphLevel, leafNodes, modes, highLevelModes, ...
                                                            options, fileList)
    display('Vocabulary learning has started.');                          
    %% ========== Step 0: Set initial data structures ==========
    vocabulary = cell(options.maxLevels,1);
    mainGraph = cell(options.maxLevels,1);
    allOppositeModes = cell(options.maxLevels,1);
    
    %% ========== Step 1: Create first vocabulary and graph layers with existing node/edge info ==========
    %% Step 1.1: Prepare intermediate data structures for sequential processing.
    vocabulary(1) = {vocabLevel};
    mainGraph(1) = {graphLevel};
    previousModes = [];
    
    %% Get number of valid images in which we get gabor responses.
    numberOfImages = numel(unique([graphLevel.imageId]));
    
    %% Print first vocabulary and graph level.
    if options.debug
        visualizeLevel( vocabulary{1}, 1, previousModes, 0, options);
        visualizeImages( fileList, graphLevel, leafNodes, 1, options, 'train' );
    end
    
    %% Calculate statistics from this graph.
    display('........ Estimating statistics for level 1..');
    [avgShareability, avgCoverage] = saveStats(vocabLevel, graphLevel, leafNodes, numberOfImages, options, 'preInhibition', 1);
    display(['........ Average coverage of leaf nodes: ' num2str(avgCoverage) ', while average shareability is: ' num2str(avgShareability) '.']); 
    
    %% ========== Step 2: Infer new parts by discovering frequent subs in data. ==========
    for levelItr = 2:options.maxLevels
        %% Step 2.0: Get opposite edge information, if there is any.
        oppositeModes = [];
        if ~isempty(modes) && numel(modes)>=(levelItr-1)
            currentModes = modes{levelItr-1};
            numberOfModes = size(currentModes,1);
            oppositeModes = zeros(numberOfModes,1);
            for modeItr = 1:numberOfModes
               oppositeMode = currentModes(modeItr,:);
               oppositeMode(1:2) = oppositeMode(2:-1:1);
               oppositeMode(3:4) = oppositeMode(3:4) * -1;
               [~, oppositeModeIdx] = ismember(oppositeMode, currentModes, 'rows' );
               oppositeModes(modeItr) = oppositeModeIdx;
            end
        else
            currentModes = [];
        end
        allOppositeModes(levelItr-1) = {oppositeModes};
        
        %% Step 2.1: Run knowledge discovery to learn frequent compositions.
        [vocabLevel, graphLevel] = discoverSubs(vocabLevel, graphLevel, oppositeModes, ...
            options, options.currentFolder, false, levelItr-1);
        
        %% If no new subs have been found, finish processing.
        if isempty(vocabLevel)
           % Write previous level's appearances to the output folder.
           vocabulary = vocabulary(1:(levelItr-1),:);
           mainGraph = mainGraph(1:(levelItr-1),:);
           allOppositeModes = allOppositeModes(1:(levelItr-1),:);
           break; 
        end
        
        %% Step 2.2: If necessary, find realizations in graph level again. 
        % This is necessary since 'exe' type implementation does not always 
        % calculate mdl correctly.  
        if strcmp(options.subdue.implementation, 'exe')
            % Reorder substructres so that the ones having higher mdl score are
            % actually on top. The mdl score here is not really MDL, but rather
            % a metric depending on size * frequency.
            [vocabLevel, graphLevel] = reorderSubs(vocabLevel, graphLevel);

            %% Assign realizations R of next graph level (l+1), and fill in their bookkeeping info.
            previousLevel = mainGraph{levelItr-1};
            graphLevel = fillBasicInfo(previousLevel, graphLevel, leafNodes);

            %% Inhibition! We process the current level to eliminate some of the nodes in the final graph.
            % The rules here are explained in the paper. Basically, each node
            % should introduce a novelty (cover a new area of the image). If it
            % fails to do so, the one which has a lower mdl ranking is
            % discarded. *Natural selection*
            
            [graphLevel] = applyLocalInhibition(vocabLevel, graphLevel, currentModes, options, levelItr);
            [remainingComps, ~, ~] = unique([graphLevel.labelId], 'stable');

            % Eliminate unused compositions from vocabulary.
            vocabLevel = vocabLevel(1, remainingComps);

            %% To get more accurate realizations, search words in vocabulary one by one, again.
            % Then, we will apply the local inhibition function again.
            graphLevel = collectInstances(vocabLevel, mainGraph{levelItr-1}, oppositeModes, options, levelItr-1);

            % Fill in basic info, again.
            graphLevel = fillBasicInfo(previousLevel, graphLevel, leafNodes);  
        else
            %% Assign realizations R of next graph level (l+1), and fill in their bookkeeping info.
            previousLevel = mainGraph{levelItr-1};
            graphLevel = fillBasicInfo(previousLevel, graphLevel, leafNodes);
        end
        
        %% Calculate statistics from this graph.
        display('........ Before we apply inhibition, estimating statistics..');
        [avgShareability, avgCoverage] = saveStats(vocabLevel, graphLevel, leafNodes, numberOfImages, options, 'preInhibition', levelItr);
        display(['........ Average Coverage: ' num2str(avgCoverage) ', average shareability of compositions: ' num2str(avgShareability) '.']); 
        
        %% Inhibition! We process the current level to eliminate some of the nodes in the final graph.
        % The rules here are explained in the paper. Basically, each node
        % should introduce a novelty (cover a new area of the image). If it
        % fails to do so, the one which has a lower mdl ranking is
        % discarded. *Natural selection*
        display('........ Applying inhibition.');
        [graphLevel] = applyLocalInhibition(vocabLevel, graphLevel, currentModes, options, levelItr);
        [remainingComps, ~, IC] = unique([graphLevel.labelId], 'stable');

        % Learn number of elements, so that we can order the vocabulary
        % based on number of occurences.
        
        
        % Assign new labels of the remaining realizations.
        IC = num2cell(IC);
        [graphLevel.labelId] = deal(IC{:});

        % Eliminate unused compositions from vocabulary.
        vocabLevel = vocabLevel(1, remainingComps);
        
        % Reorder substructres so that the ones having higher mdl score are
        % actually on top. The mdl score here is not really MDL, but rather
        % a metric depending on size * frequency.
%        [vocabLevel, graphLevel] = reorderSubs(vocabLevel, graphLevel);
        
        %% Calculate statistics from this graph.
        [avgShareability, avgCoverage] = saveStats(vocabLevel, graphLevel, leafNodes, numberOfImages, options, 'postInhibition', levelItr);
        
        % display debugging info.
        display(['........ Inhibition applied with novelty thr: ' num2str(options.noveltyThr) ' and edge novelty thr: ' num2str(options.edgeNoveltyThr) '.']);
        display(['........ Remaining: ' num2str(numel(graphLevel)) ' realizations belonging to ' num2str(numel(vocabLevel)) ' compositions.']);
        display(['........ Average Coverage: ' num2str(avgCoverage) ', average shareability of compositions: ' num2str(avgShareability) '.']); 
        
        % Set the sign of all nodes to 1. When negative graphs are introduced,
        % this part should CHANGE.
        [graphLevel.sign] = deal(1);
        
        %% Step 2.4: Create the parent relationships between current level and previous level.
        vocabulary = mergeIntoGraph(vocabulary, vocabLevel, leafNodes, levelItr, 0);
        mainGraph = mergeIntoGraph(mainGraph, graphLevel, leafNodes, levelItr, 1);
        
        %% Step 2.5: Create object graphs G_(l+1) for the next level, l+1.
        % Extract the edges between new realizations to form the new object graphs.
        [modes, highLevelModes, mainGraph] = extractEdges(mainGraph, options, levelItr, modes, highLevelModes);
        graphLevel = mainGraph{levelItr};
        
        %% Print vocabulary and graph level to output images (reconstruction).
        if options.debug
           display('........ Visualizing previous level...');
           if ~isempty(vocabLevel)
               visualizeLevel( vocabLevel, levelItr, modes{levelItr-1}, numel(vocabulary{levelItr-1}), options);
               visualizeImages( fileList, graphLevel, leafNodes, levelItr, options, 'train' );
           end
        end
        
        %% Step 2.6: If no new edges found, kill program.
        newEdgesAvailable = ~isempty(cat(1, mainGraph{levelItr}.adjInfo));
        if ~newEdgesAvailable
            vocabulary = vocabulary(1:(levelItr),:);
            mainGraph = mainGraph(1:(levelItr),:);
            allOppositeModes = allOppositeModes(1:(levelItr-1),:);
        end
    end
end