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
function [nextVocabLevel, nextGraphLevel, optimalThreshold, isSupervisedSelectionRunning, previousAccuracy] = runSubdue(vocabLevel, ...
    graphLevel, orgThreshold, nodeDistanceMatrix, edgeDistanceMatrix, categoryArrIdx, validationIdx, ...
    supervisedSelectionFlag, isSupervisedSelectionRunning, previousAccuracy, level1Nodes, levelItr, options)
    %% First thing we do is to convert vocabLevel and graphLevel into different data structures.
    % This process is done to assure fast, vectorized operations.
    % Initialize the priority queue.
    bestSubs = [];
    parentSubs = [];
    nextVocabLevel = [];
    nextGraphLevel = [];
    
    %% Get the parameters.
    evalMetric = options.subdue.evalMetric;
    mdlNodeWeight = options.subdue.mdlNodeWeight;
    mdlEdgeWeight = options.subdue.mdlEdgeWeight;
    isMDLExact = options.subdue.isMDLExact;
    overlap = options.subdue.overlap;
    beam = options.subdue.beam;
    dynamicBeam = Inf;
    nsubs = options.subdue.nsubs;
    maxTime = options.subdue.maxTime;
    minSize = options.subdue.minSize;
    maxSize = options.subdue.maxSize;
    isSupervised = options.subdue.supervised;
    maxThreshold = options.subdue.maxThreshold;
    minThreshold = options.subdue.minThreshold;
    singlePrecision = options.singlePrecision;
    optimizationFlag = options.optimizationFlag;
    partSelectionFlag = options.partSelectionFlag;
    optimalThreshold = orgThreshold;
    minRFCoverage = options.missingNodeThr;
    if options.validationFlag
        validationFolds = options.validationFolds;
    else
        validationFolds = 1; 
    end
    if optimizationFlag
        singleNodeThreshold = maxThreshold + singlePrecision;
    else
        singleNodeThreshold = orgThreshold +singlePrecision; % Hard threshold for cost of matching two subs.
    end
    parentsPerSet = 50;
    
    % At this point we get more subs than we need, since we're trying to
    % optimize based on the number of subs.
    numberOfReconstructiveSubs = options.reconstruction.numberOfReconstructiveSubs;
    
    %% Initialize data structures.
    display('[SUBDUE] Initializing data structures for internal use..');
    % Helper data structures.
    assignedEdges = {graphLevel.adjInfo};    
    if isempty(assignedEdges) 
        return;
    else
        allEdgeNodePairs = cat(1,assignedEdges{:});
        if isempty(allEdgeNodePairs)
            return;
        end
        numberOfAllEdges = size(allEdgeNodePairs,1);
        clear allEdgeNodePairs;
    end
    nonemptyEdgeIdx = cellfun(@(x) ~isempty(x), assignedEdges);
    allLabels = cat(1, graphLevel.labelId);
    assignedEdges(nonemptyEdgeIdx) = cellfun(@(x) [x(:,1:3), allLabels(x(:,2))], ...
        assignedEdges(nonemptyEdgeIdx), 'UniformOutput', false);
    allEdges(numel(graphLevel)) = AdjInfo();
    [allEdges.adjInfo] = deal(assignedEdges{:});
    allEdgeNodePairs = cat(1,allEdges.adjInfo);
    allEdgeNodePairs = allEdgeNodePairs(:,1:2);
    clear assignedEdges;
    
    % Get leaf nodes for each child.
    allLeafNodes = {graphLevel.leafNodes};
    
    % Obtain real node labels and edge distributions for likelihood
    % processing.
    realNodeLabels = [graphLevel.realLabelId];
    nodePositions = cat(1, graphLevel.position);
    edgeCoords = options.edgeCoords;
    
     % Get number of probabilistic choices for vocabulary nodes.
    numberOfProbabilisticChoicesArr = [vocabLevel.numberOfProbabilisticChoices];
    [~, validVocabLevel, ~] = unique([vocabLevel.label], 'stable');
    numberOfProbabilisticChoicesArr = numberOfProbabilisticChoicesArr(validVocabLevel);
    
    % If no edges are present, time to return.
    allSigns = uint8(cat(1, graphLevel.sign));
    imageIdx = cat(1, graphLevel.imageId);
    categoryArrIdx = uint8(categoryArrIdx(cat(1, graphLevel.imageId)))';
    validationIdx = validationIdx(cat(1, graphLevel.imageId));
    
    % Learn possible leaf nodes within every RF, and save them for future
    % use.
    possibleLeafNodes = cell(numel(allLeafNodes),1);
%     bounds = calculateRFBounds(nodePositions, levelItr, options, false);
%     for nodeItr = 1:size(allSigns,1);
%         possibleLeafNodes{nodeItr} = int32(find(level1Nodes(:,1) == imageIdx(nodeItr) & ...
%              level1Nodes(:,2) >= bounds(nodeItr,1) & level1Nodes(:,2) <= bounds(nodeItr,3) & ...
%              level1Nodes(:,3) >= bounds(nodeItr,2) & level1Nodes(:,3) <= bounds(nodeItr,4)));
%     end
     for nodeItr = 1:size(allSigns,1)
          if isempty(allEdges(nodeItr).adjInfo)
               possibleLeafNodes{nodeItr} = (graphLevel(nodeItr).leafNodes)';
          else
               possibleLeafNodes{nodeItr} = unique(cat(2, graphLevel(unique(allEdges(nodeItr).adjInfo(:,1:2))).leafNodes))';
          end
     end
    
%     % Check script.
%     for nodeItr = 1:size(allSigns,1);
%           adjInfo=allEdges(nodeItr).adjInfo;
%           if isempty(adjInfo)
%                checkNodes = graphLevel(nodeItr).leafNodes';
%           else
%                checkNodes = unique(cat(2, graphLevel(unique(allEdges(nodeItr).adjInfo(:,1:2))).leafNodes))';
%           end
%           if nnz(ismember(checkNodes, possibleLeafNodes{nodeItr})) ~= numel(checkNodes)
%                nodeItr
%           end
%      end
                  
    % Graph size formulation.
    if strcmp(evalMetric, 'likelihood')
         graphSize = numberOfAllEdges + numel(graphLevel);
    else
         % Graph size formulation is very simple: edgeWeight * #edges + edgeWeight * #nodes. 
         graphSize = numberOfAllEdges * mdlEdgeWeight + ...
             numel(graphLevel) * mdlNodeWeight;
    end
    
    % Find the total cost of matching considering the max size, and set the
    % threshold.
   if optimizationFlag
        adaptiveThreshold = maxThreshold * (single(maxSize) * 2 - 1) + singlePrecision;
   else
        adaptiveThreshold = orgThreshold * (single(maxSize) * 2 - 1) + singlePrecision;
   end
   
    %% Step 1:Find single-vertex subs, and put them into beamSubs.
    display('[SUBDUE] Creating single node subs..');
    singleNodeSubs = getSingleNodeSubs(allLabels, allSigns, ...
        nodeDistanceMatrix, categoryArrIdx, validationIdx, adaptiveThreshold);
    parentSubs = addToQueue(singleNodeSubs, parentSubs, numel(singleNodeSubs));
    singleNodeSubsFinal = [];
    
    %% Step 2: Main loop
    startTime = tic;
    % Let's find the adaptive threshold with the size of the minimum 
    currentSize = 2;
    if optimizationFlag
        overallThreshold = maxThreshold * (single(maxSize) * 2 - 1) + singlePrecision;
    else
        overallThreshold = orgThreshold * (single(maxSize) * 2 - 1) + singlePrecision;
    end
    minMdlScoreFinal = -inf; % This two will increase as more and more subs are discovered.
    bestMdlScoresFinal = [];
    
    while ~isempty(parentSubs)
        minMdlScoreExtend = -inf;   
        % This two will increase as more and more subs are discovered.
        % If we're building a reconstructive hierarchy, we do not want to
        % eliminate subs, as we may need them later. That's why we are
        % setting adaptiveMinThreshold to the minimum value.
        % adaptiveThreshold is used to eliminate the subs to be put in the
        % final node list.
        if optimizationFlag
            adaptiveThreshold = maxThreshold * (single(currentSize) * 2 - 1) + singlePrecision;
            adaptiveMinThreshold = minThreshold * (single(currentSize) * 2 - 1) + singlePrecision;
        else
            adaptiveThreshold = orgThreshold * (single(currentSize) * 2 - 1) + singlePrecision;
            adaptiveMinThreshold = adaptiveThreshold;
        end
        
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
                 singleNodeSubsFinal = getSingleNodeSubs(allLabels, allSigns, ...
                    nodeDistanceMatrix, categoryArrIdx, validationIdx, singleNodeThreshold);
               
               % Evaluate them.
                singleNodeSubsFinal = evaluateSubs(singleNodeSubsFinal, [], evalMetric, allEdges, allEdgeNodePairs, ...
                    allSigns, graphSize, overlap, mdlNodeWeight, mdlEdgeWeight, isMDLExact, isSupervised, false, ...
                    allLeafNodes, minRFCoverage, possibleLeafNodes);

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
            parfor parentItr = processedSet
                %% Step 2.2: Extend head in all possible directions into childSubs.
                display(['[SUBDUE/Parallel] Expanding sub ' num2str(parentItr) ' of size ' num2str(currentSize-1) '..']);
                childSubs = extendSub(parentSubs(parentItr), allEdges, nodeDistanceMatrix, edgeDistanceMatrix, ...
                     overallThreshold);
                if isempty(childSubs) 
                    continue;
                end

                %% Step 2.3: Eliminate subs that match to more frequent subs.
                 [childSubs, ~] = getNonOverlappingSubs(childSubs, nodeDistanceMatrix, edgeDistanceMatrix, ...
                      adaptiveMinThreshold, singlePrecision);

                 % Get the list of final subs and indices of subs chosen
                 % for extension. 
                [childSubsFinal, childSubsExtend] = getFinalSubs(childSubs, adaptiveThreshold);
                 
                %% Step 2.4: Evaluate childSubs, find their instances.
                [childSubsFinal, childSubsExtend, validSubs, validExtSubs] = evaluateSubs(childSubsFinal, childSubsExtend, evalMetric, allEdges, allEdgeNodePairs, ...
                    allSigns, graphSize, overlap, mdlNodeWeight, mdlEdgeWeight, isMDLExact, isSupervised, ...
                    false, allLeafNodes, minRFCoverage, possibleLeafNodes);
                
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
            
            %% Based on the best mdl scores, update minMdlScoreFinal and minMdlScoreExtend.
            newMdlScores = [bestMdlScoresFinal, [mdlScoreArrFinal{:}]];
            if numel(newMdlScores) > nsubs
                newMdlScores = sort(newMdlScores, 'descend');
                minMdlScoreFinal = newMdlScores(nsubs);
            end
            newMdlScores = [mdlScoreArrExtend{:}];
            if numel(newMdlScores) > dynamicBeam
                newMdlScores = sort(newMdlScores, 'descend');
                minMdlScoreExtend = newMdlScores(dynamicBeam);
            end
            
            % In addition, remove the subs that have low mdl scores.
            % First, we handle final subs.
            if minMdlScoreFinal ~= -inf
                validBestSubIdx = bestMdlScoresFinal >= minMdlScoreFinal;
                bestSubs = bestSubs(validBestSubIdx);
                nonemptyArrIdx = cellfun(@(x) ~isempty(x), childSubArrFinal);
                childSubArrFinal(nonemptyArrIdx) = cellfun(@(x, y) x(y >= minMdlScoreFinal), childSubArrFinal(nonemptyArrIdx), mdlScoreArrFinal(nonemptyArrIdx), 'UniformOutput', false);
                mdlScoreArrFinal(nonemptyArrIdx) = cellfun(@(x) x(x >= minMdlScoreFinal), mdlScoreArrFinal(nonemptyArrIdx), 'UniformOutput', false);
            end
            
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
        if ~isempty(childSubArrFinal)
            %% A substructure has many different parse trees at this point.
            % We remove duplicate subs from childSubArr.
            [childSubArrFinal] = removeDuplicateSubs(childSubArrFinal);
            
            % Check for minSize to put to final sub array.
            if currentSize >= minSize
                bestSubs = addToQueue(childSubArrFinal, bestSubs, nsubs);
                if numel(bestSubs) == nsubs
                    minMdlScoreFinal = min([bestSubs.mdlScore]);
                end
            end
        end
        if ~isempty(childSubArrExtend)
            %% A substructure has many different parse trees at this point.
            % We remove duplicate subs from childSubArr.
            [childSubArrExtend] = removeDuplicateSubs(childSubArrExtend);
            
            % In this step, we take into account where each node has been
            % seen. We don't really want to eliminate subs that encompass
            % children seen only once.
            % Sort subs.
            mdlScoresTemp = [childSubArrExtend.mdlScore];
            [~, idx] = sort(mdlScoresTemp, 'descend');
            childSubArrExtend = childSubArrExtend(idx);
            
            % Write down children-s parents and only keep one parent for
            % each child.
            encounterArr = inf(numel(graphLevel), 1);
            leafNodeCountArr = zeros(numel(graphLevel),1);
            for itr = numel(childSubArrExtend):-1:1
                 instanceCenterIdx = childSubArrExtend(itr).instanceCenterIdx;
                 instanceChildren = childSubArrExtend(itr).instanceChildren;
                 
                 %% We employ a greedy approach, and keep info of instances with maximum number of leaf nodes.
                 coveredLeafNodes = cell(numel(instanceCenterIdx),1);
                 for childItr = 1:numel(instanceCenterIdx)
                      coveredLeafNodes{childItr} = int32(unique(cat(2, allLeafNodes{instanceChildren(childItr,:)})));
                 end
                 coveredLeafNodeCount = cellfun(@(x) numel(x), coveredLeafNodes);
                 
                 % Greedy elimination. Check current number of covered leaf
                 % nodes for a center. If this is more, update the relevant
                 % sub indices.
                 comparisonNodeCounts = leafNodeCountArr(instanceCenterIdx);
                 assgnIdx = coveredLeafNodeCount > comparisonNodeCounts; 
                 leafNodeCountArr(instanceCenterIdx(assgnIdx)) = coveredLeafNodeCount(assgnIdx);
                 encounterArr(instanceCenterIdx(assgnIdx)) = itr;
                 
                 % Save sub id.
                 encounterArr(instanceCenterIdx) = min(itr, encounterArr(instanceCenterIdx));
            end
            childSubArrExtend = childSubArrExtend(setdiff(unique(encounterArr), Inf));
            
            % Check for maxSize to put to extendedSubs (parentSubs for next level).
            if currentSize < maxSize
                extendedSubs = addToQueue(childSubArrExtend, extendedSubs, beam);
            end
        end
            
        clear childSubArrExtend childSubArrFinal;
        
        %% Step 2.6: Swap expandedSubs with parentSubs.
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
        bestSubs = addToQueue(singleNodeSubsFinal, bestSubs, nsubs + numel(singleNodeSubsFinal)); 
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
    allRemainingChildren = unique(cat(1, allRemainingChildren{:}));
    maxCoverage = numel(allRemainingChildren) / numel(graphLevel);
    display(['[SUBDUE] Maximal coverage possible: %' num2str(maxCoverage * 100) '.']);
    clear vocabLevel graphLevel
    
  %% If required, we'll pick best parts based on the reconstruction of the data.
    if numberOfInstances>0
        display(['[SUBDUE] We have found ' num2str(numberOfBestSubs) ' subs with ' num2str(numberOfInstances) ' instances.']);
       if partSelectionFlag
           [selectedSubs, selectedThreshold, optimalAccuracy] = selectParts(bestSubs, realNodeLabels,...
               nodePositions, edgeCoords, singlePrecision, ...
               numberOfReconstructiveSubs, orgThreshold, ...
               validationFolds, validationIdx, categoryArrIdx, imageIdx, allSigns, allLeafNodes, possibleLeafNodes, ...
               isSupervisedSelectionRunning);
           
           % If supervision flag is set, and the performance has dropped
           % since the previous iteration, we switch to supervision.
           if supervisedSelectionFlag && ~isSupervisedSelectionRunning && ...
               optimalAccuracy < previousAccuracy
           
               display(['[SUBDUE] Validation accuracy has dropped from %' num2str(100 * previousAccuracy) ' to %' num2str(100 * optimalAccuracy) '.']);
               display('[SUBDUE] We switch to supervised learning from now on.');
               isSupervisedSelectionRunning = true;
               [selectedSubs, selectedThreshold, previousAccuracy] = selectParts(bestSubs, realNodeLabels, ...
                   nodePositions, edgeCoords, singlePrecision, ...
                   numberOfReconstructiveSubs, orgThreshold, ...
                   validationFolds, validationIdx, categoryArrIdx, imageIdx, allSigns, allLeafNodes, possibleLeafNodes, ...
                   isSupervisedSelectionRunning);
           else
               previousAccuracy = optimalAccuracy;
           end
           bestSubs = selectedSubs;
           optimalThreshold = selectedThreshold;
           
           % Re-evaluate best subs.
           bestSubs = evaluateSubs(bestSubs, [], 'mdl', allEdges, allEdgeNodePairs, ...
               allSigns, graphSize, overlap, mdlNodeWeight, mdlEdgeWeight, 1, isSupervised, false, ...
               allLeafNodes, minRFCoverage, possibleLeafNodes);
           
           % Sort bestSubs by their mdl scores.
           mdlScores = [bestSubs.mdlScore];
           [~, mdlSortIdx] = sort(mdlScores, 'descend');
           bestSubs = bestSubs(mdlSortIdx);
       end
       
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
        instanceChildrenMappings = cell(numberOfBestSubs,1);
        parfor bestSubItr = 1:numberOfBestSubs
           instanceChildren = bestSubs(bestSubItr).instanceChildren;
           instanceMappings = bestSubs(bestSubItr).instanceMappings;
           tempNum = size(instanceChildren,1);
           
           % Find children descriptors and save 'em.
           childrenDescriptors = zeros(tempNum, maxSubSize, 'int32');
           childrenMappings = zeros(tempNum, maxSubSize, 'uint8');
           childrenDescriptors(:, 1:size(instanceChildren,2)) = instanceChildren;
           childrenMappings(:, 1:size(instanceChildren,2)) = instanceMappings;
           instanceChildrenDescriptors{bestSubItr} = childrenDescriptors;
           instanceChildrenMappings{bestSubItr} = childrenMappings;
        end
        instanceChildrenDescriptors = cat(1, instanceChildrenDescriptors{:});
        instanceChildrenMappings = cat(1, instanceChildrenMappings{:});

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
           vocabLevel(bestSubItr).normMdlScore = bestSubs(remainingBestSubs(bestSubItr)).normMdlScore;

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
           
           % Assign the number of choices up to this point. 
           vocabLevel(bestSubItr).numberOfProbabilisticChoices = sum(numberOfProbabilisticChoicesArr(vocabLevel(bestSubItr).children)) + (2 * int32(numberOfEdges) + 1); %#ok<PFBNS>
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
               
               % Get mappings and signs for assignment.
                instanceMappings = instanceChildrenMappings(actualInstanceOffset:(actualInstanceOffset + numberOfInstances-1),:);
                instanceMappings( :, all(~instanceMappings,1) ) = [];
                instanceChildren = mat2cell(instanceChildren, ones(numberOfInstances,1), size(instanceChildren,2));
                instanceMappings = mat2cell(instanceMappings, ones(numberOfInstances,1), size(instanceMappings,2));
                instanceSigns = num2cell(allSigns(centerIdx));
               
               % Calculate activations.
               instanceActivations = ones(numberOfInstances,1, 'single'); 
               instanceActivations = num2cell(instanceActivations);
               
               % Assign fields to graphs.
               [graphLevel(actualInstanceOffset:(actualInstanceOffset + numberOfInstances-1)).labelId] = deal(labelIds{:});
               [graphLevel(actualInstanceOffset:(actualInstanceOffset + numberOfInstances-1)).realLabelId] = deal(labelIds{:});
               [graphLevel(actualInstanceOffset:(actualInstanceOffset + numberOfInstances-1)).children] = deal(instanceChildren{:});
               [graphLevel(actualInstanceOffset:(actualInstanceOffset + numberOfInstances-1)).mapping] = deal(instanceMappings{:});
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
    
   clear allEdges allEdgeNodePairs allSigns bestSubs remainingInstanceLabels instanceChildrenDescriptors;
   nextVocabLevel = vocabLevel;
   nextGraphLevel = graphLevel;
   clear vocabLevel graphLevel;
end