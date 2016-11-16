%> Name: runSubdue
%>
%> Description: This function runs SUBDUE (self) with the input image and
%> parameters defined in options. 
%>
%> @param vocabLevel Input vocabulary level. If preDefinedSearch is 0,
%> ignored. If 1, compositions in this vocabulary level are detected in
%> graphLevel.
%> @param graphLevel The current object graphs' level. The graphLevel's
%> nodes should be sorted by their image id. 
%> label1 encodes exact opposite geometric information of each mode in its 
%> index, and same applies to every row of this list. If empty, simply ignored.
%> @param options Program options.
%>
%> Explanation of the algorithm:
%> 
%> Given Graph V = (N, E), with node list N and edge Eist E, 
%>       Queue = [], ExtendQueue = [], FinalQueue = [], options.
%> Queue, ExtendQueue, FinalQueue are ordered by mdl scores. They are
%> automatically trimmed to options.beam elements after any addition.
%> 1- Find single node subs in N and add them to Queue. Occurences of each
%> node type is counted, and Queue is ordered by the number of occurences.
%> 
%> 2- For each sub S in Queue,
%>      2.1- Extend S by exploring all possible extensions in G 
%>      into ChildSubs.
%> 
%>      2.2- Remove any duplicate subs from ChildSubs. 
%> 
%>      2.3- Evauate each sub in ChildSubs and assign their mdl scores.
%>
%>      2.4- Add each sub in ChildSubs into ExtendQueue.
%>
%> 3- Add each sub in ExtendQueue to FinalQueue. No duplicate subs are
%> added to FinalQueue if they already exist.
%> 
%> 4- Write final subs and their instances to vocabLevel, graphLevel.
%>
%> @retval nextVocabLevel Next vocabulary level ([] if pre-defined search
%> is on). 
%> @retval nextGraphLevel Next object graphs' level.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 05.02.2014
%> Ver 1.1 on 01.09.2014 Removal of global parameters.
%> Ver 1.2 on 02.09.2014 Adding display commentary.
function [nextVocabLevel, nextGraphLevel, isSupervisedSelectionRunning, previousAccuracy] = runSubdue(vocabLevel, ...
    graphLevel, level1Coords, categoryArrIdx, supervisedSelectionFlag, isSupervisedSelectionRunning, previousAccuracy, levelItr, options)
    %% First thing we do is to convert vocabLevel and graphLevel into different data structures.
    % This process is done to assure fast, vectorized operations.
    % Initialize the priority queue.
    bestSubs = [];
    parentSubs = [];
    nextVocabLevel = [];
    nextGraphLevel = [];

    if levelItr < 3
        singleInstanceFlag = false; 
        minSize = 2;
    else
        minSize = options.subdue.minSize;
        singleInstanceFlag = true;
    end
     %    numberOfPrevSubs = numel(vocabLevel);

     % Obtain initial coverage.
     initialCoverage = numel(fastsortedunique(sort(cat(2, graphLevel.leafNodes))));
    
    %% Get the parameters.
    evalMetric = options.subdue.evalMetric;
    mdlNodeWeight = options.subdue.mdlNodeWeight;
    mdlEdgeWeight = options.subdue.mdlEdgeWeight;
    isMDLExact = options.subdue.isMDLExact;
    overlap = options.subdue.overlap;
    beam = options.subdue.beam;
    nsubs = options.subdue.nsubs;
    maxTime = options.subdue.maxTime;
    maxSize = options.subdue.maxSize;
    singleNodeSubThreshold = options.subdue.singleNodeSubThreshold;
    maxShareability = options.maxShareability;
    halfRFSize = ceil(options.receptiveFieldSize/2);
    smallHalfMatrixSize = (options.smallReceptiveFieldSize+1)/2;
    if ismember(levelItr, options.smallRFLayers)
        rfRadius = smallHalfMatrixSize;
    else
        rfRadius = halfRFSize;
    end
    
    % If this is category level, we change the minimum RF coverage.
    if levelItr + 1 == options.categoryLevel
        minRFCoverage = options.categoryLevelCoverage;
    else
        minRFCoverage = options.missingNodeThr;
    end
    parentsPerSet = 400;
    
%     % On lower levels, we do not need to extend too much.
%     maxSize = min(maxSize, levelItr + 1);
    
    % At this point we get more subs than we need, since we're trying to
    % optimize based on the number of subs.
    numberOfReconstructiveSubs = options.reconstruction.numberOfReconstructiveSubs;
    
    %% Create folder for visualizations.
    folderName = [options.debugFolder '/level' num2str(levelItr+1) '/discovery'];
    if ~exist(folderName, 'dir')
         mkdir(folderName);
    end
    
    %% Initialize data structures.
    display('[SUBDUE] Initializing data structures for internal use..');
    % Helper data structures.
    allEdges = {graphLevel.adjInfo};    
    if isempty(allEdges) 
        return;
    else
        allEdgeNodePairs = cat(1,allEdges{:});
        if isempty(allEdgeNodePairs)
            return;
        end
        clear allEdgeNodePairs;
    end
    nonemptyEdgeIdx = cellfun(@(x) ~isempty(x), allEdges);
    allLabels = cat(1, graphLevel.labelId);
    allEdges(nonemptyEdgeIdx) = cellfun(@(x) [x(:,1:3), allLabels(x(:,2))], ...
        allEdges(nonemptyEdgeIdx), 'UniformOutput', false);
    allEdgeNodePairs = cat(1,allEdges{:});
    allEdgeNodePairs = allEdgeNodePairs(:,1:2);
    allEdgeCounts = hist(double(allEdgeNodePairs(:,1)), 1:numel(graphLevel))';
    allCoords = cat(1, graphLevel.position);
    clear assignedEdges;
    avgDegree = size(allEdgeNodePairs,1)/numel(graphLevel);
    
    % Get leaf nodes for each child.
    allLeafNodes = {graphLevel.leafNodes};
    
    % If no edges are present, time to return.
    allSigns = uint8(cat(1, graphLevel.sign));
    imageIdx = cat(1, graphLevel.imageId);
    
%     % Build a very large (sparse) possible adjacency matrix.
%     display('[SUBDUE] Creating linkage for nodes to be used in finding cliques.');
%     maxImageId = max(imageIdx);
%     tempArr = cell(maxImageId,1);
%     for imageItr = 1:maxImageId
%          nodeListIdx = find(imageIdx == imageItr);
%          if nnz(nodeListIdx) == 0
%               continue;
%          end
%          
%          % If there are nodes, look for nearby ones.
%          distances = pdist(allCoords(nodeListIdx, :)) < rfRadius;
%          for nodeItr = 1:numel(nodeListIdx)
%               adjNodes = find(distances(nodeItr,:));
%               
%          end
%     end
    
    adjMatrix = sparse(double(allEdgeNodePairs(:,1)), double(allEdgeNodePairs(:,2)), ones(size(allEdgeNodePairs,1),1) > 0, numel(graphLevel), numel(graphLevel));
    
    
    % Learn possible leaf nodes within every RF, and save them for future
    % use. For this one, we consider only the nodes within the receptive
    % field.
    % Learn stride.
    if strcmp(options.filterType, 'gabor')
         stride = options.gabor.stride;
    else
         stride = options.auto.stride;
    end
    poolDim = options.poolDim;
    
    % Calculate pool factor.
    if levelItr> 1
        poolFactor = nnz(~ismembc(2:levelItr, options.noPoolingLayers));
    else
        poolFactor = 0;
    end
    
    % Allocate space for counts for every center, and find number of
    % accessible leaf nodes for each receptive field.
    level1CoordsPooled = calculatePooledPositions(level1Coords(:,2:3), poolFactor, poolDim, stride);
    possibleLeafNodeCounts = zeros(numel(allLeafNodes),1);
    if minRFCoverage > 0
         parfor nodeItr = 1:size(allSigns,1)
             % Obtain node related data.
             nodeEdges = allEdges{nodeItr};
             nodeCoords = allCoords(nodeItr,:);
             nodeChildren = allLeafNodes{nodeItr};
             
             % Calculate number of accessible leaf nodes within the RF.
             if isempty(nodeEdges)
                  possibleLeafNodeCounts(nodeItr)  = numel(nodeChildren);
             else
                 tempChildren = cat(1, nodeItr, nodeEdges(:,2));
                 tempLeafNodes = fastsortedunique(sort(cat(2, allLeafNodes{tempChildren})));
                 tempLeafNodeCoords = level1CoordsPooled(tempLeafNodes,:);
                 
                 % Keep only leaf nodes which exist within this RF.
                 tempLeafNodes = tempLeafNodes(tempLeafNodeCoords(:,1) > nodeCoords(1) - rfRadius & ...
                     tempLeafNodeCoords(:,1) < nodeCoords(1) + rfRadius & ...
                     tempLeafNodeCoords(:,2) > nodeCoords(2) - rfRadius & ...
                     tempLeafNodeCoords(:,2) < nodeCoords(2) + rfRadius);
                 
                 % Save the number.
                 possibleLeafNodeCounts(nodeItr) = numel(tempLeafNodes);
             end
         end
    end
   
    %% Step 1:Find single-vertex subs, and put them into beamSubs.
    display('[SUBDUE] Creating single node subs..');
    singleNodeSubs = getSingleNodeSubs(allLabels, allSigns);
    parentSubs = addToQueue(singleNodeSubs, parentSubs, numel(singleNodeSubs));
    singleNodeSubsFinal = [];
    
    %% Step 2: Main loop
    startTime = tic;
    currentSize = 2;
    bestMdlScoresFinal = [];
    
    while ~isempty(parentSubs)
        minMdlScoreFinal = -inf; % This two will increase as more and more subs are discovered.
        minMdlScoreExtend = -inf;   
        if ~isempty(bestSubs)
            bestMdlScoresFinal = [bestSubs.mdlScore];
        end
        
        % Check time. If it took too long, end processing. 
        % Check parent's size. If it is too large, end processing.
        elapsedTime = toc(startTime);
        if size(parentSubs(1).edges,1) >= (maxSize-1)
            display(['[SUBDUE] Maximum size of ' num2str(maxSize) ' nodes reached, terminating search.']);
            break;
        end
        if elapsedTime > maxTime 
            display(['[SUBDUE] Time is up! Breaking from process just before extending subs of size ' num2str(size(parentSubs(1).edges,1)+1) '..']);
            break;
        end
        
        %% Evaluate single node subs, if required.
        if minSize == 1
            if currentSize == 2
                 % Find single node subs.
                 singleNodeSubsFinal = getSingleNodeSubs(allLabels, allSigns);
               
               % Evaluate them.
                singleNodeSubsFinal = evaluateSubs(singleNodeSubsFinal, evalMetric, allEdgeCounts, allEdgeNodePairs, ...
                    allSigns, allCoords, overlap, mdlNodeWeight, mdlEdgeWeight, isMDLExact, ...
                    allLeafNodes, level1CoordsPooled, rfRadius, minRFCoverage, maxShareability, possibleLeafNodeCounts, avgDegree, singleNodeSubThreshold);

                %% Remove those with no instances. 
                centerIdxArr = {singleNodeSubsFinal.instanceCenterIdx};
                validSingleSubs = cellfun(@(x) ~isempty(x), centerIdxArr);
                singleNodeSubsFinal = singleNodeSubsFinal(validSingleSubs);

                %% Sort singleNodeSubsFinal by their scores.
                % We do not eliminate nodes with negative mdl scores,
                % as single node subs do not provide any compression.
                mdlScores = [singleNodeSubsFinal.mdlScore];
                [~, sortIdx] = sort(mdlScores, 'descend');
                singleNodeSubsFinal = singleNodeSubsFinal(sortIdx);
            end
        else
            singleNodeSubsFinal = [];
        end
        
        %% All good, execution. First, we put parents into fixed size subsets for easy parallelization.
        display(['[SUBDUE/Parallel] Expanding subs of size ' num2str(size(parentSubs(1).edges,1)+1) '..']);
        childSubArrFinal = cell(numel(parentSubs),1);
        childSubArrExtend = cell(numel(parentSubs),1);
        mdlScoreArrFinal = cell(numel(parentSubs),1);
        mdlScoreArrExtend = cell(numel(parentSubs),1);
        numberOfParentSubs = numel(parentSubs);
        if numberOfParentSubs > parentsPerSet
            setDistributions = 1:parentsPerSet:numberOfParentSubs;
            if setDistributions(end) ~= numberOfParentSubs
                setDistributions = [setDistributions, numberOfParentSubs]; %#ok<AGROW>
            end
            setDistributions = setDistributions(2:(end)) - setDistributions(1:(end-1));
            setDistributions(end) = setDistributions(end) + 1;
        else
            setDistributions = numberOfParentSubs;
        end
        parentSubSets = cell(mat2cell(1:numel(parentSubs), 1, setDistributions));
        preservedSubs = [];
        for setItr = 1:numel(parentSubSets)
            %% Step 2.1: If it has been too long, we need to finish execution.
            elapsedTime = toc(startTime);
            haltedParent = -1;
            if elapsedTime > maxTime && currentSize > 2
                haltedParent = parentSubSets{setItr}(1);
                display(['[SUBDUE] Time is up! Breaking from process before ' ...
                    'extending parent ' num2str(haltedParent) ' of size ' num2str(size(parentSubs(1).edges,1)+1) '..']);
                break;
            end
            
            %% All good, continue with the main algorithm.
            processedSet = parentSubSets{setItr};
            display(['[SUBDUE/Parallel] Expanding set ' num2str(setItr) '/' num2str(numel(parentSubSets)) ' of size ' num2str(currentSize-1) ', containing ' num2str(numel(processedSet)) ' subs..']);
            try
                 parfor parentItr = processedSet
                     %% Step 2.2: Extend head in all possible directions into childSubs.
      %               display(['[SUBDUE/Parallel] Expanding sub ' num2str(parentItr) ' of size ' num2str(currentSize-1) '..']);
                     childSubs = extendSub(parentSubs(parentItr), allEdges, allEdgeCounts, singleInstanceFlag);
                     if isempty(childSubs) 
                         continue;
                     end

                      % Get the list of final subs and indices of subs chosen
                      % for extension. 
                     childSubsFinal = childSubs;
                     childSubsExtend = childSubs;

                     %% Step 2.4: Evaluate childSubs, find their instances.
                     [childSubsFinal, validSubs, validExtSubs] = evaluateSubs(childSubsFinal, evalMetric, allEdgeCounts, allEdgeNodePairs, ...
                         allSigns, allCoords, overlap, mdlNodeWeight, mdlEdgeWeight, isMDLExact, ...
                         allLeafNodes, level1CoordsPooled, rfRadius, minRFCoverage, maxShareability, possibleLeafNodeCounts, avgDegree, singleNodeSubThreshold);

                     % Assign mdl scores of subs chosen for extension as well. 
                     [childSubsFinal, childSubsExtend] = copyMdlScores(childSubsFinal, childSubsExtend);
                     childSubsFinal = childSubsFinal(validSubs);
                     childSubsExtend = childSubsExtend(validExtSubs);

                     %% Eliminate childSubs which will render useless to preserve memory.
                     if ~isempty(childSubsFinal)
                          finalMdlScores = [childSubsFinal.mdlScore];
                          validMdlScoreIdxFinal = finalMdlScores > minMdlScoreFinal;
                          childSubsFinal = childSubsFinal(validMdlScoreIdxFinal);
                     end
                     if ~isempty(childSubsExtend)
                          mdlScores = [childSubsExtend.mdlScore];
                          validMdlScoreIdxExtend = mdlScores > minMdlScoreExtend;
                          childSubsExtend = childSubsExtend(validMdlScoreIdxExtend);
                     end

                     %% Sort childSubs by the mdl scores.
                     if ~isempty(childSubsFinal)
                          mdlScores = [childSubsFinal.mdlScore];
                          [sortedMdlScoresFinal, sortIdx] = sort(mdlScores, 'descend');
                          childSubsFinal = childSubsFinal(sortIdx);
                     else
                          sortedMdlScoresFinal = [];
                     end
                     if ~isempty(childSubsExtend)
                          mdlScores = [childSubsExtend.mdlScore];
                          [sortedMdlScoresExtend, sortIdx] = sort(mdlScores, 'descend');
                          childSubsExtend = childSubsExtend(sortIdx);
                     else
                          sortedMdlScoresExtend = [];
                     end

                     %% Save childSubs and extended subs.
                     childSubArrFinal(parentItr) = {childSubsFinal};
                     childSubArrExtend(parentItr) = {childSubsExtend};
                     mdlScoreArrFinal(parentItr) = {sortedMdlScoresFinal};
                     mdlScoreArrExtend(parentItr) = {sortedMdlScoresExtend};
                 end
            catch
                 for parentItr = processedSet
                     %% Step 2.2: Extend head in all possible directions into childSubs.
      %               display(['[SUBDUE/Parallel] Expanding sub ' num2str(parentItr) ' of size ' num2str(currentSize-1) '..']);
                     childSubs = extendSub(parentSubs(parentItr), allEdges, allEdgeCounts, singleInstanceFlag);
                     if isempty(childSubs) 
                         continue;
                     end

                      % Get the list of final subs and indices of subs chosen
                      % for extension. 
                     childSubsFinal = childSubs;
                     childSubsExtend = childSubs;

                     %% Step 2.4: Evaluate childSubs, find their instances.
                     [childSubsFinal, validSubs, validExtSubs] = evaluateSubs(childSubsFinal, evalMetric, allEdgeCounts, allEdgeNodePairs, ...
                         allSigns, allCoords, overlap, mdlNodeWeight, mdlEdgeWeight, isMDLExact, ...
                         allLeafNodes, level1CoordsPooled, rfRadius, minRFCoverage, maxShareability, possibleLeafNodeCounts, avgDegree, singleNodeSubThreshold);

                     % Assign mdl scores of subs chosen for extension as well. 
                     [childSubsFinal, childSubsExtend] = copyMdlScores(childSubsFinal, childSubsExtend);
                     childSubsFinal = childSubsFinal(validSubs);
                     childSubsExtend = childSubsExtend(validExtSubs);

                     %% Eliminate childSubs which will render useless to preserve memory.
                     if ~isempty(childSubsFinal)
                          finalMdlScores = [childSubsFinal.mdlScore];
                          validMdlScoreIdxFinal = finalMdlScores > minMdlScoreFinal;
                          childSubsFinal = childSubsFinal(validMdlScoreIdxFinal);
                     end
                     if ~isempty(childSubsExtend)
                          mdlScores = [childSubsExtend.mdlScore];
                          validMdlScoreIdxExtend = mdlScores > minMdlScoreExtend;
                          childSubsExtend = childSubsExtend(validMdlScoreIdxExtend);
                     end

                     %% Sort childSubs by the mdl scores.
                     if ~isempty(childSubsFinal)
                          mdlScores = [childSubsFinal.mdlScore];
                          [sortedMdlScoresFinal, sortIdx] = sort(mdlScores, 'descend');
                          childSubsFinal = childSubsFinal(sortIdx);
                     else
                          sortedMdlScoresFinal = [];
                     end
                     if ~isempty(childSubsExtend)
                          mdlScores = [childSubsExtend.mdlScore];
                          [sortedMdlScoresExtend, sortIdx] = sort(mdlScores, 'descend');
                          childSubsExtend = childSubsExtend(sortIdx);
                     else
                          sortedMdlScoresExtend = [];
                     end

                     %% Save childSubs and extended subs.
                     childSubArrFinal(parentItr) = {childSubsFinal};
                     childSubArrExtend(parentItr) = {childSubsExtend};
                     mdlScoreArrFinal(parentItr) = {sortedMdlScoresFinal};
                     mdlScoreArrExtend(parentItr) = {sortedMdlScoresExtend};
                 end
           end
                 
            
            %% Based on the best mdl scores, update minMdlScoreFinal and minMdlScoreExtend.
            %% TODO TO BE OPENED ON DISCRIMINATIVE LEARNING>
%             newMdlScores = [mdlScoreArrFinal{:}];
%             if numel(newMdlScores) > nsubs
%                 newMdlScores = sort(newMdlScores, 'descend');
%                 minMdlScoreFinal = newMdlScores(nsubs);
%             end
          
%            encounterArr = inf(numel(graphLevel), 1);
%            for itr = numel(childSubArrFinal):-1:1
%                 instanceCenterIdx = childSubArrFinal(itr).instanceCenterIdx;
%                 encounterArr(instanceCenterIdx) = itr;
%            end
%            keptSubs = setdiff(unique(encounterArr), Inf)';
%            childSubArrFinal = childSubArrFinal(keptSubs);

            % In addition, remove the subs that have low mdl scores.
            % First, we handle final subs.
            preservedSubs = cat(2, preservedSubs, childSubArrFinal{processedSet});
            childSubArrFinal(processedSet) = {[]};
            
            % Remove unnecessary subs.
            encounterArr = inf(numel(graphLevel), 1);
            for itr = numel(preservedSubs):-1:1
                instanceCenterIdx = preservedSubs(itr).instanceCenterIdx;
                encounterArr(instanceCenterIdx) = itr;
            end
            keptSubs = setdiff(unique(encounterArr), Inf)';
            preservedSubs = preservedSubs(keptSubs);
            
%             if minMdlScoreFinal ~= -inf
%                 nonemptyArrIdx = cellfun(@(x) ~isempty(x), childSubArrFinal);
%                 childSubArrFinal(nonemptyArrIdx) = cellfun(@(x, y) x(y >= minMdlScoreFinal), childSubArrFinal(nonemptyArrIdx), mdlScoreArrFinal(nonemptyArrIdx), 'UniformOutput', false);
%                 mdlScoreArrFinal(nonemptyArrIdx) = cellfun(@(x) x(x >= minMdlScoreFinal), mdlScoreArrFinal(nonemptyArrIdx), 'UniformOutput', false);
%             end
            
            % Then, we handle extended subs.
            if minMdlScoreExtend ~= -inf
                nonemptyArrIdx = cellfun(@(x) ~isempty(x), childSubArrExtend);
                childSubArrExtend(nonemptyArrIdx) = cellfun(@(x, y) x(y >= minMdlScoreExtend), childSubArrExtend(nonemptyArrIdx), mdlScoreArrExtend(nonemptyArrIdx), 'UniformOutput', false);
                mdlScoreArrExtend(nonemptyArrIdx) = cellfun(@(x) x(x >= minMdlScoreExtend), mdlScoreArrExtend(nonemptyArrIdx), 'UniformOutput', false);
            end
        end
        
        %% Add each children group in childGroupArr into extendedSubs.
        display('[SUBDUE] Merging all children and putting them into bestSubs if they match final criteria.');
        extendedSubs = [];
        %% Step 2.4: Add childSubs to extendedSubs and bestSubs.
        childSubArrFinal = cat(2,childSubArrFinal{:});
        childSubArrExtend = cat(2, childSubArrExtend{:});
        % Add children to both queues.
        if ~isempty(preservedSubs)
            %% A substructure has many different parse trees at this point.
            % We remove duplicate subs from childSubArr.
%            [childSubArrFinal] = removeDuplicateSubs(childSubArrFinal, numberOfPrevSubs);

           % Remove unnecessary subs.
           encounterArr = inf(numel(graphLevel), 1);
           for itr = numel(childSubArrFinal):-1:1
                instanceCenterIdx = childSubArrFinal(itr).instanceCenterIdx;
                encounterArr(instanceCenterIdx) = itr;
           end
           keptSubs = setdiff(unique(encounterArr), Inf)';
           childSubArrFinal = childSubArrFinal(keptSubs);
            
            % Remove excess subs.
            preservedSubs = addToQueue(preservedSubs, [], nsubs);
            
            % Check for minSize to put to final sub array.
            if currentSize >= minSize
                bestSubs = addToQueue(preservedSubs, bestSubs, nsubs * (currentSize - 1));
            end
        end
        
        if ~isempty(childSubArrExtend)
            %% A substructure has many different parse trees at this point.
            % We remove duplicate subs from childSubArr.
%            [childSubArrExtend] = removeDuplicateSubs(childSubArrExtend, numberOfPrevSubs);
            
            % In this step, we take into account where each node has been
            % seen. We don't really want to eliminate subs that encompass
            % children seen only once.
            % Sort subs.
            mdlScoresTemp = [childSubArrExtend.mdlScore];
            [~, idx] = sort(mdlScoresTemp, 'descend');
            childSubArrExtend = childSubArrExtend(idx);
            
            % Write down children-s parents and only keep one parent for
            % each child.
 %           if numel(childSubArrExtend) > beam
                encounterArr = inf(numel(graphLevel), 1);
                for itr = numel(childSubArrExtend):-1:1
                     instanceCenterIdx = childSubArrExtend(itr).instanceCenterIdx;
                     encounterArr(instanceCenterIdx) = itr;
                end
                childSubArrExtend = childSubArrExtend(setdiff(unique(encounterArr), Inf));
%            end
            
            % Check for maxSize to put to extendedSubs (parentSubs for next level).
            if currentSize < maxSize
                extendedSubs = addToQueue(childSubArrExtend, extendedSubs, beam);
            end
        end
            
        clear childSubArrExtend childSubArrFinal;
        
        %% Step 2.6: Swap expandedSubs with parentSubs, if we're still allowed to proceed.
        if haltedParent == -1
            display('[SUBDUE] Swapping children with parents and going on..');
            parentSubs = extendedSubs;
            currentSize = currentSize + 1;
            clear extendedSubs;
        else
            break;
        end
    end
    
    %% At this point, we add single node subs to the list of known subs.
    if ~isempty(singleNodeSubsFinal)
        bestSubs = addToQueue(singleNodeSubsFinal, bestSubs, numel(bestSubs) + numel(singleNodeSubsFinal)); 
    end
    %% Step 3: Create nextVocabLevel and nextGraphLevel from bestSubs.
    numberOfBestSubs = numel(bestSubs);
    
    % If no subs detected, return.
    if numberOfBestSubs < 1
       return;
    end
    
    % Remove duplicate instances from each sub.
    display('[SUBDUE/Parallel] Removing duplicate instances from each sub, and updating their match scores..');
    bestSubs = removeDuplicateInstances(bestSubs);
    
    %% Allocate space for new graphLevel and vocabLevel.
    numberOfInstances = 0;
    for bestSubItr = 1:numberOfBestSubs
        numberOfInstances = numberOfInstances + size(bestSubs(bestSubItr).instanceCenterIdx,1);
    end
    
    %% Calculate maximum coverage possible based on the data at hand.
    allRemainingChildren = {bestSubs.instanceChildren};
    allRemainingChildren = cellfun(@(x) x(:), allRemainingChildren, 'UniformOutput', false);
    allRemainingChildren = fastsortedunique(sort(cat(1, allRemainingChildren{:})));
    maxCoverage = numel(fastsortedunique(sort(cat(2, allLeafNodes{allRemainingChildren})))) / initialCoverage;
    display(['[SUBDUE] Maximal coverage possible: %' num2str(maxCoverage * 100) '.']);
    clear vocabLevel graphLevel
    
  %% If required, we'll pick best parts based on the reconstruction of the data.
    if numberOfInstances>0
       display(['[SUBDUE] We have found ' num2str(numberOfBestSubs) ' subs with ' num2str(numberOfInstances) ' instances.']);
       [selectedSubs, optimalAccuracy] = selectParts(bestSubs, ...
           numberOfReconstructiveSubs, categoryArrIdx, imageIdx, allSigns, allLeafNodes, level1Coords,...
           isSupervisedSelectionRunning, folderName);

       % If supervision flag is set, and the performance has dropped
       % since the previous iteration, we switch to supervision.
       if supervisedSelectionFlag && ~isSupervisedSelectionRunning && ...
           optimalAccuracy < previousAccuracy

           display(['[SUBDUE] Validation accuracy has dropped from %' num2str(100 * previousAccuracy) ' to %' num2str(100 * optimalAccuracy) '.']);
           display('[SUBDUE] We switch to supervised learning from now on.');
           isSupervisedSelectionRunning = true;
           [selectedSubs, previousAccuracy] = selectParts(bestSubs, ...
           numberOfReconstructiveSubs, categoryArrIdx, imageIdx, allSigns, allLeafNodes, isSupervisedSelectionRunning, folderName);
       else
           previousAccuracy = optimalAccuracy;
       end
      bestSubs = selectedSubs;

%        % Re-evaluate best subs.
%        bestSubs = evaluateSubs(bestSubs, 'mdl', allEdgeCounts, allEdgeNodePairs, ...
%            allSigns, allCoords, overlap, mdlNodeWeight, mdlEdgeWeight, false, ...
%            allLeafNodes, level1CoordsPooled, rfRadius, minRFCoverage, maxShareability, possibleLeafNodeCounts, avgDegree, singleNodeSubThreshold);

       % Sort bestSubs by their mdl scores.
       mdlScores = [bestSubs.mdlScore];
       [~, mdlSortIdx] = sort(mdlScores, 'descend');
       bestSubs = bestSubs(mdlSortIdx);
       
       %% Learn activations for instances, and save the children for future use.
        display('[SUBDUE/Parallel] Calculating activations for each instance, and saving instance nodes for later use..');
        numberOfBestSubs = numel(bestSubs);
        % Learn the maximum children count in instances.
        maxSubSize = 1;
        for bestSubItr = 1:numberOfBestSubs
            maxSubSize = max(maxSubSize, (size(bestSubs(bestSubItr).edges,1)+1));
        end

        % Find the children for each instance, and use it as a descriptor.
        % We'll then find unique descriptors, which is the initial phase of
        % inhibition.
        instanceChildrenDescriptors = cell(numberOfBestSubs,1);
        parfor bestSubItr = 1:numberOfBestSubs
           instanceChildren = bestSubs(bestSubItr).instanceChildren;
           tempNum = size(instanceChildren,1);
           
           % Find children descriptors and save 'em.
           childrenDescriptors = zeros(tempNum, maxSubSize, 'int32');
           childrenDescriptors(:, 1:size(instanceChildren,2)) = instanceChildren;
           instanceChildrenDescriptors{bestSubItr} = childrenDescriptors;
        end
        instanceChildrenDescriptors = cat(1, instanceChildrenDescriptors{:});

        remainingBestSubs = 1:numberOfBestSubs;
        IA = 1:size(instanceChildrenDescriptors,1);
        numberOfInstances = numel(IA);
        
        %% Fill in vocabLevel and graphLevel.
       %Allocate space for new graph/vocab level.
       display('[SUBDUE] Processing has finished. Writing all to output and quitting.');
       vocabLevel(numberOfBestSubs) = VocabNode();
       graphLevel(numberOfInstances) = GraphNode();
       
        % First, we start with vocabLevel.
        parfor bestSubItr = 1:numberOfBestSubs
           % Assign label of sub.
           vocabLevel(bestSubItr).label = int32(bestSubItr);

           % Assign mdl score.
           vocabLevel(bestSubItr).mdlScore = bestSubs(remainingBestSubs(bestSubItr)).mdlScore; %#ok<PFBNS>

           % Assign sub edges (definition)
           numberOfEdges = size(bestSubs(remainingBestSubs(bestSubItr)).edges,1);
           if numberOfEdges > 0
               subEdges = [ones(numberOfEdges,1), ...
                   (2:1:(numberOfEdges+1))', bestSubs(remainingBestSubs(bestSubItr)).edges(:,1), ones(numberOfEdges,1)];
               vocabLevel(bestSubItr).adjInfo = subEdges;
               vocabLevel(bestSubItr).children = [bestSubs(remainingBestSubs(bestSubItr)).centerId; bestSubs(remainingBestSubs(bestSubItr)).edges(:,2)]';
           else
               vocabLevel(bestSubItr).children = bestSubs(remainingBestSubs(bestSubItr)).centerId;
           end
        end
        
        % Now, we fill in the info of vocabLevel.
        actualInstanceOffset = 1;
        instanceOffset = 1;
        remainingSubItr = 1;
        for bestSubItr = 1:numel(bestSubs)
           if ~ismember(bestSubItr,remainingBestSubs)
               instanceOffset = instanceOffset + numel(bestSubs(bestSubItr).instanceCenterIdx);
           else
               %% Assign instances to this sub.
               instanceEndOffset = instanceOffset + numel(bestSubs(bestSubItr).instanceCenterIdx) - 1;
               remainingInstanceIdx = (IA(IA>= instanceOffset & IA <= instanceEndOffset) - instanceOffset) + 1;
               numberOfInstances = numel(remainingInstanceIdx);
               labelIds = num2cell(repmat(int32(remainingSubItr), numberOfInstances,1));

               %% Create required fields such as centerIdx, edges, children and sign array to assign to instances.
               centerIdx = bestSubs(bestSubItr).instanceCenterIdx(remainingInstanceIdx,:);
               % Get graph match costs.
               numberOfInstances = numel(centerIdx);
               
               % Get instance children
               instanceChildren = instanceChildrenDescriptors(actualInstanceOffset:(actualInstanceOffset + numberOfInstances-1),:);
               instanceChildren( :, all(~instanceChildren,1) ) = [];
               numberOfInstances = numel(centerIdx);
               instanceChildren = mat2cell(instanceChildren, ones(numberOfInstances,1), size(instanceChildren,2));
               instanceSigns = num2cell(allSigns(centerIdx));
               
               % Calculate activations.
               instanceActivations = ones(numberOfInstances,1, 'single'); 
               instanceActivations = num2cell(instanceActivations);
               
               % Assign fields to graphs.
               [graphLevel(actualInstanceOffset:(actualInstanceOffset + numberOfInstances-1)).labelId] = deal(labelIds{:});
               [graphLevel(actualInstanceOffset:(actualInstanceOffset + numberOfInstances-1)).realLabelId] = deal(labelIds{:});
               [graphLevel(actualInstanceOffset:(actualInstanceOffset + numberOfInstances-1)).children] = deal(instanceChildren{:});
               [graphLevel(actualInstanceOffset:(actualInstanceOffset + numberOfInstances-1)).sign] = deal(instanceSigns{:});
               [graphLevel(actualInstanceOffset:(actualInstanceOffset + numberOfInstances-1)).activation] = deal(instanceActivations{:});
               instanceOffset = instanceEndOffset + 1;
               actualInstanceOffset = actualInstanceOffset + numberOfInstances;
               remainingSubItr = remainingSubItr + 1;
           end
        end
    else
        vocabLevel = [];
        graphLevel = [];
    end
    
   nextVocabLevel = vocabLevel;
   nextGraphLevel = graphLevel;
   clearvars -except nextVocabLevel nextGraphLevel isSupervisedSelectionRunning previousAccuracy
end