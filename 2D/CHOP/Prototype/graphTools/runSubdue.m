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
function [nextVocabLevel, nextGraphLevel, optimalThreshold] = runSubdue(vocabLevel, graphLevel, nodeDistanceMatrix, edgeDistanceMatrix, categoryArrIdx, options)
    %% First thing we do is to convert vocabLevel and graphLevel into different data structures.
    % This process is done to assure fast, vectorized operations.
    % Initialize the priority queue.
    
    bestSubs = [];
    parentSubs = [];
    nextVocabLevel = [];
    nextGraphLevel = [];
    allLeafNodes = {graphLevel.leafNodes};
    prevGraphNodes = unique([graphLevel.leafNodes]);
    prevGraphNodeCount = numel(prevGraphNodes);
    maxPrevGraphNodeId = max(prevGraphNodes);
    
    %% Get the parameters.
    evalMetric = options.subdue.evalMetric;
    mdlNodeWeight = options.subdue.mdlNodeWeight;
    mdlEdgeWeight = options.subdue.mdlEdgeWeight;
    isMDLExact = options.subdue.isMDLExact;
    overlap = options.subdue.overlap;
    beam = options.subdue.beam;
    nsubs = options.subdue.nsubs;
    maxTime = options.subdue.maxTime;
    minSize = options.subdue.minSize;
    maxSize = options.subdue.maxSize;
    isSupervised = options.subdue.supervised;
    orgThreshold = single(options.subdue.threshold);
    maxThreshold = options.subdue.maxThreshold;
    minThreshold = options.subdue.minThreshold;
    maxDepth = options.subdue.thresholdSearchMaxDepth;
    singlePrecision = options.singlePrecision;
    reconstructionFlag = options.reconstruction.flag;
    stoppingCoverage = options.reconstruction.stoppingCoverage;
    optimalThreshold = orgThreshold;
    if reconstructionFlag
        singleNodeThreshold = maxThreshold + singlePrecision;
    else
        singleNodeThreshold = orgThreshold +singlePrecision; % Hard threshold for cost of matching two subs.
    end
    
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
    
    % If no edges are present, time to return.
    allSigns = uint8(cat(1, graphLevel.sign));
    categoryArrIdx = uint8(categoryArrIdx(cat(1, graphLevel.imageId)))';
    
    % Graph size formulation is very simple: edgeWeight * #edges + edgeWeight * #nodes. 
    graphSize = numberOfAllEdges * mdlEdgeWeight + ...
        numel(graphLevel) * mdlNodeWeight;
    
    %% Step 1:Find single-vertex subs, and put them into beamSubs.
    display('[SUBDUE] Creating single node subs..');
    singleNodeSubs = getSingleNodeSubs(allLabels, allSigns, ...
        nodeDistanceMatrix, categoryArrIdx, singleNodeThreshold);
%    parentSubs = addToQueue(singleNodeSubs, parentSubs, options.subdue.beam);
    parentSubs = addToQueue(singleNodeSubs, parentSubs, numel(singleNodeSubs));
    singleNodeSubsFinal = [];
    
    %% Step 2: Main loop
    startTime = tic;
    % Let's find the adaptive threshold with the size of the minimum 
    currentSize = 2;
    minMdlScore = -inf; % This two will increase as more and more subs are discovered.
    bestMdlScores = [];
    
    while ~isempty(parentSubs)
        % If we're building a reconstructive hierarchy, we do not want to
        % eliminate subs, as we may need them later. That's why we are
        % setting adaptiveMinThreshold to the minimum value.
        if reconstructionFlag
            adaptiveThreshold = maxThreshold * (single(currentSize) * 2 - 1) + singlePrecision;
            adaptiveMinThreshold = minThreshold * (single(currentSize) * 2 - 1) + singlePrecision;
 %           adaptiveMinThreshold = minThreshold + singlePrecision;
        else
            adaptiveThreshold = orgThreshold * (single(currentSize) * 2 - 1) + singlePrecision;
            adaptiveMinThreshold = adaptiveThreshold;
%            adaptiveMinThreshold = orgThreshold + singlePrecision;
        end
        
        if ~isempty(bestSubs)
            bestMdlScores = [bestSubs.mdlScore];
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
                singleNodeSubsFinal = evaluateSubs(parentSubs, evalMetric, allEdges, allEdgeNodePairs, ...
                    allSigns, graphSize, overlap, mdlNodeWeight, mdlEdgeWeight, isMDLExact, isSupervised);

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
        
        %% All good, execution.
        display(['[SUBDUE/Parallel] Expanding subs of size ' num2str(size(parentSubs(1).edges,1)+1) '..']);
        childSubArr = cell(numel(parentSubs),1);
        mdlScoreArr = cell(numel(parentSubs),1);
        numberOfParentSubs = numel(parentSubs);
        if numberOfParentSubs > 50
            setDistributions = 1:50:numberOfParentSubs;
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
            if elapsedTime > maxTime
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
                childSubs = extendSub(parentSubs(parentItr), allEdges, nodeDistanceMatrix, edgeDistanceMatrix, 1, adaptiveThreshold);
                if isempty(childSubs) 
                    continue;
                end
                
                %% Step 2.3: Eliminate subs that match to more frequent subs.
                 [childSubs, ~] = getNonOverlappingSubs(childSubs, nodeDistanceMatrix, edgeDistanceMatrix, ...
                     adaptiveMinThreshold, singlePrecision);

                %% Step 2.4: Evaluate childSubs, find their instances.
                % Setting overlap to 1, since we've already checked for it
                % in extendSub.
                childSubs = evaluateSubs(childSubs, evalMetric, allEdges, allEdgeNodePairs, ...
                    allSigns, graphSize, overlap, mdlNodeWeight, mdlEdgeWeight, isMDLExact, isSupervised);
                
                %% Eliminate childSubs which will render useless to preserve memory.
                mdlScores = [childSubs.mdlScore];
                validMdlScoreIdx = mdlScores > minMdlScore;
                childSubs = childSubs(validMdlScoreIdx);
                    
                % If no subs are remaining, continue.
                if isempty(childSubs) 
                    continue;
                end
                
                %% Sort childSubs and childSubs by the mdl scores.
                mdlScores = [childSubs.mdlScore];
                [sortedMdlScores, sortIdx] = sort(mdlScores, 'descend');
                childSubs = childSubs(sortIdx);
                
                %% Save childSubs and extended subs.
                childSubArr(parentItr) = {childSubs};
                mdlScoreArr(parentItr) = {sortedMdlScores};
            end
            %% Based on the best mdl scores, update minMdlScore.
            newMdlScores = [bestMdlScores, [mdlScoreArr{:}]];
            if numel(newMdlScores) > nsubs
                newMdlScores = sort(newMdlScores, 'descend');
                minMdlScore = newMdlScores(nsubs);
            end
            
            % In addition, remove the subs that have low mdl scores.
            if minMdlScore ~= -inf
                validBestSubIdx = bestMdlScores >= minMdlScore;
                bestSubs = bestSubs(validBestSubIdx);
                nonemptyArrIdx = cellfun(@(x) ~isempty(x), childSubArr);
                childSubArr(nonemptyArrIdx) = cellfun(@(x, y) x(y >= minMdlScore), childSubArr(nonemptyArrIdx), mdlScoreArr(nonemptyArrIdx), 'UniformOutput', false);
                mdlScoreArr(nonemptyArrIdx) = cellfun(@(x) x(x >= minMdlScore), mdlScoreArr(nonemptyArrIdx), 'UniformOutput', false);
            end
        end
        
        %% Add each children group in childGroupArr into extendedSubs.
        display('[SUBDUE] Merging all children and putting them into bestSubs if they match final criteria.');
        extendedSubs = [];
        %% Step 2.4: Add childSubs to extendedSubs and bestSubs.
        childSubArr = cat(2,childSubArr{:});
        % Add children to the both queues.
        if ~isempty(childSubArr)
            %% A substructure has many different parse trees at this point.
            % We find matching parse trees, and combine their instances.
            numberOfChildSubs = numel(childSubArr);
            [childSubArr, childSubArrToEvaluate] = removeDuplicateSubs(childSubArr);
            
            % If the set of subs to be re-evaluated is not empty, we
            % evaluate them again, and 
            if ~isempty(childSubArrToEvaluate)
                childSubArrToEvaluate = evaluateSubs(childSubArrToEvaluate, evalMetric, allEdges, allEdgeNodePairs, ...
                    allSigns, graphSize, overlap, mdlNodeWeight, mdlEdgeWeight, isMDLExact, isSupervised);
                childSubArr = addToQueue(childSubArrToEvaluate, childSubArr, numberOfChildSubs);
            end
            
            % Check for maxSize to put to extendedSubs (parentSubs for next level).
            if currentSize < maxSize
                extendedSubs = addToQueue(childSubArr, extendedSubs, beam);
            end
            % Check for minSize to put to final sub array.
            if currentSize >= minSize
                bestSubs = addToQueue(childSubArr, bestSubs, nsubs);
                if numel(bestSubs) == nsubs
                    minMdlScore = min([bestSubs.mdlScore]);
                end
            end
        end
        clear childSubArr childSubArrToEvaluate;
        
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
    
    display('[SUBDUE] Processing has finished. Writing all to output and quitting.');
    %% Step 3: Create nextVocabLevel and nextGraphLevel from bestSubs.
    numberOfBestSubs = numel(bestSubs);
    
    % If no subs detected, return.
    if numberOfBestSubs < 1
       return;
    end
    
    %% Allocate space for new graphLevel and vocabLevel.
    numberOfInstances = 0;
    for bestSubItr = 1:numberOfBestSubs
        numberOfInstances = numberOfInstances + size(bestSubs(bestSubItr).instanceCenterIdx,1);
    end
    clear vocabLevel graphLevel
    if numberOfInstances>0
       %% If required, we'll pick best parts based on the reconstruction of the data.
       if reconstructionFlag
           [bestSubs, optimalThreshold] = getReconstructiveParts(bestSubs, ...
               allLeafNodes, prevGraphNodeCount, maxPrevGraphNodeId, nodeDistanceMatrix, edgeDistanceMatrix, singlePrecision, ...
               stoppingCoverage, numberOfReconstructiveSubs, minThreshold, maxThreshold, maxDepth);
           bestSubs = evaluateSubs(bestSubs, 'mdl', allEdges, allEdgeNodePairs, ...
                    allSigns, graphSize, 1, mdlNodeWeight, mdlEdgeWeight, isMDLExact, isSupervised); 
           for bestSubItr = 1:numel(bestSubs)
               bestSubs(bestSubItr).mdlScore = bestSubs(bestSubItr).mdlScore / numel(bestSubs(bestSubItr).instanceCenterIdx);
           end
           mdlScores = [bestSubs.mdlScore];
           [~, mdlSortIdx] = sort(mdlScores, 'descend');
           bestSubs = bestSubs(mdlSortIdx);
       end
       
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
        instanceMatchCosts = cell(numberOfBestSubs,1);
        remainingInstanceLabels = cell(numberOfBestSubs,1);
        parfor bestSubItr = 1:numberOfBestSubs
           instanceChildren = bestSubs(bestSubItr).instanceChildren;
           numberOfInstances = size(instanceChildren,1);
           % Find children descriptors and save 'em.
           childrenDescriptors = zeros(numberOfInstances, maxSubSize, 'int32');
           childrenDescriptors(:, 1:size(instanceChildren,2)) = instanceChildren;
           instanceChildrenDescriptors{bestSubItr} = childrenDescriptors;
           instanceMatchCosts{bestSubItr} = bestSubs(bestSubItr).instanceMatchCosts;
           
           % Find instance labels.
           instanceLabels = repmat(int32(bestSubItr), numberOfInstances, 1);
           remainingInstanceLabels{bestSubItr} = instanceLabels;
        end
        instanceChildrenDescriptors = cat(1, instanceChildrenDescriptors{:});
        instanceMatchCosts = cat(1, instanceMatchCosts{:});
        remainingInstanceLabels = cat(1, remainingInstanceLabels{:});
        [~, sortIdx] = sort(instanceMatchCosts, 'ascend');
        orderedInstanceChildrenDescriptors = instanceChildrenDescriptors(sortIdx,:);

       %% Update data structures by selecting unique children.
       % Find unique rows, which correspond to unique instances.
       [~, IA, ~] = unique(orderedInstanceChildrenDescriptors, 'rows', 'stable');
       IA = sort(sortIdx(IA));
       
       % Get the remaining sub count.
       remainingBestSubs = unique(remainingInstanceLabels(IA))';
       numberOfBestSubs = numel(remainingBestSubs);
       numberOfInstances = numel(IA);
       instanceChildrenDescriptors = instanceChildrenDescriptors(IA,:);
       
       %Allocate space for new graph/vocab level.
       vocabLevel(numberOfBestSubs) = options.vocabNode;
       graphLevel(numberOfInstances) = options.graphNode;

        %% Fill in vocabLevel and graphLevel.
        % First, we start with vocabLevel.
        for bestSubItr = 1:numberOfBestSubs
           % Assign label of sub.
           vocabLevel(bestSubItr).label = int32(bestSubItr);

           % Assign mdl score.
           vocabLevel(bestSubItr).mdlScore = bestSubs(remainingBestSubs(bestSubItr)).mdlScore;
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
               numberOfInstances = numel(centerIdx);
               instanceChildren = instanceChildrenDescriptors(actualInstanceOffset:(actualInstanceOffset + numberOfInstances-1),:);
               instanceChildren( :, all(~instanceChildren,1) ) = [];
               instanceChildren = mat2cell(instanceChildren, ones(numberOfInstances,1), size(instanceChildren,2));
               instanceSigns = num2cell(allSigns(centerIdx));

               % Assign fields to graphs.
               [graphLevel(actualInstanceOffset:(actualInstanceOffset + numberOfInstances-1)).labelId] = deal(labelIds{:});
               [graphLevel(actualInstanceOffset:(actualInstanceOffset + numberOfInstances-1)).children] = deal(instanceChildren{:});
               [graphLevel(actualInstanceOffset:(actualInstanceOffset + numberOfInstances-1)).sign] = deal(instanceSigns{:});
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

%> Name: getSingleNodeSubs
%>
%> Description: getSingleNodeSubs is used to obtain single-node subs from
%> labels of all instances in allLabels. The result is a number of
%> subs representing compositions each having their instances.
%> 
%> @param allLabels Labels for every graph node.
%> @param allSigns Signs for every graph node.
%> @param categoryArrIdx Categories for every graph node.
%>
%> @retval singleNodeSubs Substructure list of single-node subs.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 24.02.2014
%> Ver 1.1 on 01.09.2014 Removal of global parameters.
function singleNodeSubs = getSingleNodeSubs(allLabels, allSigns, nodeDistanceMatrix, categoryArrIdx, threshold)
    numberOfSubs = max(allLabels);
    singleNodeSubs(numberOfSubs) = Substructure();
    validSubs = ones(numberOfSubs,1)>0;
    
    %% For each center node label type, we create a substructure.
    for subItr = 1:numberOfSubs
        distances = nodeDistanceMatrix(allLabels, subItr);
        subCenterIdx = distances < threshold;
        instances = int32(find(subCenterIdx));
        numberOfInstances = numel(instances);
        
        %Assign center id.
        singleNodeSubs(subItr).centerId = subItr;
        
        % Give maximum score so that it is at the top of the queue.
        singleNodeSubs(subItr).mdlScore = numberOfInstances;
        if numberOfInstances>0
            % Fill in instance information. 
            categoryIdx = categoryArrIdx(subCenterIdx);
            if size(categoryIdx,1) == 1
                categoryIdx = categoryIdx';
            end
            instanceSigns = allSigns(subCenterIdx,1);
            instanceDistances = distances(subCenterIdx);
            singleNodeSubs(subItr).instanceCenterIdx = instances;
            singleNodeSubs(subItr).instanceChildren = instances;
            singleNodeSubs(subItr).instanceCategories = categoryIdx;
            singleNodeSubs(subItr).instanceMatchCosts = instanceDistances;
            singleNodeSubs(subItr).instanceSigns = instanceSigns;
        else
            validSubs(subItr) = 0;
        end
    end
    
    % Eliminate those with no instances.
    singleNodeSubs = singleNodeSubs(validSubs);
end

%> Name: extendSub
%>
%> Description: extendSub(..) extends 'sub' in all possible ways
%> by extending its instances, and returns the new sub list 'extendedSubs'
%> along with instances of each sub in returned list.
%> 
%> @param sub Sub to be extended.
%> @param allEdges List of all edges in the graph (Nx4).
%>
%> @retval extendedSubs Extended sub list.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 24.02.2014
%> Ver 1.1 on 01.09.2014 Removal of global parameters.
function [extendedSubs] = extendSub(sub, allEdges, nodeDistanceMatrix, edgeDistanceMatrix, overlap, threshold)
    centerIdx = sub.instanceCenterIdx;
    subAllEdges = {allEdges(centerIdx).adjInfo}';
    % Get unused edges from sub's instances.
    allUsedEdgeIdx = sub.instanceEdges;
    if ~isempty(allUsedEdgeIdx)
        allUsedEdgeIdx = mat2cell(allUsedEdgeIdx, ones(size(allUsedEdgeIdx,1),1), size(allUsedEdgeIdx,2));
        allUnusedEdgeIdx = cellfun(@(x,y) setdiff(1:size(x,1),y), subAllEdges, allUsedEdgeIdx, 'UniformOutput', false);
        allUnusedEdges = cellfun(@(x,y) x(y,:), subAllEdges, allUnusedEdgeIdx, 'UniformOutput', false);
%        allUsedEdges = cellfun(@(x,y) x(y,:), subAllEdges, allUsedEdgeIdx, 'UniformOutput', false);
        allUnusedEdgeIdx = [allUnusedEdgeIdx{:}]';
    else
        allUnusedEdges = subAllEdges;
        allUnusedEdgeIdx = [];
%        allUsedEdges = [];
    end
    
%     % Here, we mark used nodes for this sub.
%     if ~overlap
%         validNodeArr = ones(numel(allEdges),1, 'uint8')>0;
%         validNodeArr(centerIdx) = 0;
%         if ~isempty(allUsedEdges)
%             allUsedEdges = cat(1, allUsedEdges{:});
%             allUsedEdges = allUsedEdges(:,2);
%             validNodeArr(allUsedEdges) = 0;
%         end
%     end
    
    % Record which edge belongs to which instance. 
    allEdgeInstanceIds = zeros(sum(cellfun(@(x) size(x,1), allUnusedEdges)),1);
    allEdgePrevCosts = zeros(size(allEdgeInstanceIds), 'single');
    itrOffset = 1;
    unusedEdgeCounts = cellfun(@(x) size(x,1), allUnusedEdges);
    for itr = 1:numel(allUnusedEdges)
        beginOffset = itrOffset;
        endOffset = (beginOffset+(unusedEdgeCounts(itr)-1));
        allEdgeInstanceIds(beginOffset:endOffset) = itr;
        allEdgePrevCosts(beginOffset:endOffset) = sub.instanceMatchCosts(itr);
        itrOffset = itrOffset + unusedEdgeCounts(itr);
    end         
    allUnusedEdges = cat(1, allUnusedEdges{:});
    
    % If no edges remain, exit.
    if isempty(allUnusedEdges)
        extendedSubs = [];
        return; 
    end
    
    % Remove edges leading to duplicate node use.
%     if ~overlap
%         validEdgeIdx = validNodeArr(allUnusedEdges(:,2));
%         allUnusedEdges = allUnusedEdges(validEdgeIdx,:);
%         allEdgePrevCosts = allEdgePrevCosts(validEdgeIdx);
%     end
    
    % Get unique rows of [edgeLabel, secondVertexLabel]
    uniqueEdgeTypes = unique(allUnusedEdges(:, 3:4), 'rows');
    
    % Discard any edge types already existing in sub definition.
    if ~isempty(sub.edges)
        uniqueEdgeTypes = setdiff(uniqueEdgeTypes, sub.edges, 'rows');
    end
    
    % Extend the definition in sub with each edge type in uniqueEdgeTypes.
    % In addition, we pick suitable instances, add this edge, and mark used
    % field of relevant instances.
    numberOfEdgeTypes = size(uniqueEdgeTypes,1);
    if numberOfEdgeTypes == 0
        extendedSubs = [];
        return;
    end
    
    % Save local indices for all edges, to be used later.
    allLocalEdgeIdx = cellfun(@(x) 1:size(x,1), subAllEdges, 'UniformOutput', false);
    allLocalEdgeIdx = [allLocalEdgeIdx{:}]';
%     if ~overlap
%         if size(validEdgeIdx,1) == size(allLocalEdgeIdx,1)
%             allLocalEdgeIdx = allLocalEdgeIdx(validEdgeIdx);
%         end
%         allEdgeInstanceIds = allEdgeInstanceIds(validEdgeIdx);
%         if ~isempty(allUnusedEdgeIdx)
%             allUnusedEdgeIdx = allUnusedEdgeIdx(validEdgeIdx);
%         end
%     end
    
    % Allocate space for new subs and fill them in.
    extendedSubs(numberOfEdgeTypes) = Substructure;
    for edgeTypeItr = 1:numberOfEdgeTypes
        % Assign sub-definition type and other info
        newSub = Substructure;
        newSub.edges = sub.edges;
        newSub.centerId = sub.centerId;
        newSub.mdlScore = 0;
        newSub.edges = [newSub.edges; uniqueEdgeTypes(edgeTypeItr,:)];
        
        %% Find instances of this new sub, and mark the new edges as 'used' in each of its subs.
        % Extend the subs by taking the threshold into account. Unless the
        % instance's collective matching cost surpasses the threshold, it
        % is a valid instance.
        edgesToExtendCosts = allEdgePrevCosts + ...
            edgeDistanceMatrix(allUnusedEdges(:,3), uniqueEdgeTypes(edgeTypeItr,1)) + ...
            nodeDistanceMatrix(allUnusedEdges(:,4), uniqueEdgeTypes(edgeTypeItr,2));
        edgesToExtendIdx = edgesToExtendCosts < threshold;
        
        % Remove instances leading to duplicate use.
        if ~overlap
            % We check for overlaps among center nodes/added nodes.
            idx = find(edgesToExtendIdx);
            bothNodes = allUnusedEdges(edgesToExtendIdx,1:2);
            edgesToExtendIdx(edgesToExtendIdx) = 0;
            [~, uniqueCenterNodeIdx, ~] = unique(bothNodes(:,1), 'stable');
            [~, uniqueSecondNodeIdx, ~] = unique(bothNodes(uniqueCenterNodeIdx,2), 'stable');
            idx = idx(uniqueCenterNodeIdx(uniqueSecondNodeIdx));
            edgesToExtendIdx(idx) = 1;
        end
      
        % Save instance ids.
        edgeInstanceIds = allEdgeInstanceIds(edgesToExtendIdx);
        allChildren = [sub.instanceChildren(edgeInstanceIds,:), ...
             allUnusedEdges(edgesToExtendIdx,2)];
        
        % If overlaps are not wanted, we eliminate overlapping instances
        % here. However, this elimination is not very aggressive. An
        % aggressive elimination scheme has been tried, but it tried to be
        % to strict. Useful subs have been eliminated this way. So, we have
        % decided to keep the first instance a child has been seen,
        % regardless of overlaps relating to other nodes.
%         if ~overlap
%             [~, uniqueNodeIdx ] = unique(allChildren, 'stable');
%             allChildrenValidRowsIdx = zeros(size(allChildren), 'uint8') > 0;
%             allChildrenValidRowsIdx(uniqueNodeIdx) = 1;
%             allChildrenValidRowsIdx = sum(allChildrenValidRowsIdx, 2) > 0;
%             edgeInstanceIds = edgeInstanceIds(allChildrenValidRowsIdx);
%             allChildren = allChildren(allChildrenValidRowsIdx, :);
%             idx = find(edgesToExtendIdx);
%             edgesToExtendIdx(edgesToExtendIdx) = 0;
%             idx = idx(allChildrenValidRowsIdx);
%             edgesToExtendIdx(idx) = 1;
%         end
        
        % Assign fields related to the new instances.
        newSub.instanceCenterIdx = sub.instanceCenterIdx(edgeInstanceIds);
        newSub.instanceChildren = sortrows(allChildren);
        newSub.instanceSigns = sub.instanceSigns(edgeInstanceIds);
        newSub.instanceCategories = sub.instanceCategories(edgeInstanceIds);
        newSub.instanceMatchCosts = edgesToExtendCosts(edgesToExtendIdx);
        existingEdges = sub.instanceEdges;
        combinedEdges = zeros(numel(edgeInstanceIds), size(existingEdges,2)+1, 'uint8');
        if ~isempty(existingEdges)
            existingEdges = existingEdges(edgeInstanceIds,:);
            combinedEdges(:, 1:size(existingEdges,2)) =  existingEdges;
        end
        
        % Adding the local edges to the instance edge lists.
        if ~isempty(sub.edges)
            relevantLocalEdgeIdx = allUnusedEdgeIdx(edgesToExtendIdx);
        else
            relevantLocalEdgeIdx = allLocalEdgeIdx(edgesToExtendIdx);
        end
        combinedEdges(:,end) = relevantLocalEdgeIdx;
        newSub.instanceEdges = combinedEdges;
        
        % All instances assigned, good to go.
        extendedSubs(edgeTypeItr) = newSub;
    end
end

%> Name: evaluateSubs
%>
%> Description: The evaluation function of SUBDUE. Based on the
%> evaluation metric, the value of each substructure in subs list is
%> calculated, and saved within subs. The MDL calculation takes place here.
%> Key points:
%>  1) DL estimation not rigorous. It is considered that 
%>        each node needs a label (int) and a pointer to its edges (int)
%>        each edge needs a node label (int) for its destination node, an
%>        edge label(int) and a binary isDirected label (bit)
%>     in the final graph description. Node and edge weights in DL
%>     calculation is stored in options.
%> 
%> @param subs Sub list which will be evaluated.
%> @param evalMetric Evaluation metric.
%> @param allEdges List of all edges in the graph.
%> @param allEdgeNodePairs List of all edge node pairs in the graph.
%> @param allSigns List of all signs of the nodes in the graph.
%> @param graphSize Size of the graph.
%> @param overlap If true, overlapping instances are considered in
%> evaluation of a sub.
%> @param mdlNodeWeight Node weight in DL calculations (MDL).
%> @param mdlEdgeWeight Edge weight in DL calculations (MDL).
%> @param isMDLExact If true, exact MDL calculation. Approximate otherwise.
%>
%> @retval subs Evaluated sub list.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 24.02.2014
%> Ver 1.1 on 01.09.2014 Removal of global parameters.
function [subs] = evaluateSubs(subs, evalMetric, allEdges, allEdgeNodePairs, allSigns, graphSize, overlap, mdlNodeWeight, mdlEdgeWeight, isMDLExact, isSupervised)
    numberOfSubs = numel(subs);
    parfor subItr = 1:numberOfSubs
        % Find the weight of this node, by taking the max of the category distribution. 
        if isSupervised
            categoryArr = double(subs(subItr).instanceCategories);
            weight = nnz(categoryArr == mode(categoryArr)) / numel(categoryArr);
        else
            weight = 1;
        end
        
        % We compress the object graph using the children, and the
        % edges they are involved. 
        [subScore, sub] = getSubScore(subs(subItr), allEdges, allEdgeNodePairs, evalMetric, ...
           allSigns, mdlNodeWeight, mdlEdgeWeight, ....
            overlap, isMDLExact);
        subScore = subScore * weight;
        subs(subItr) = sub;

        %% Assign the score of the sub, as well as its normalized mdl score if applicable.
        subs(subItr).mdlScore = subScore;
        if strcmp(evalMetric, 'mdl')
            subs(subItr).normMdlScore = 1 - (subScore / graphSize);
        end
    end
end

%> Name: addToQueue
%>
%> Description: Adds subs in beamQueue, which is a priority queue that orders 
%> the substructures inside by their mdlScore.
%> 
%> @param subs An array of substructures to add to beamQueue.
%> @param queue The priority queue of substructures, ordered by
%> mdlScore. Best sub has highest score.
%> @param maxSize Maximum number of elements allowed in queue.
%>
%> @retval beamQueue The priority queue of substructures, ordered by
%> mdlScore.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 06.02.2014
function [queue] = addToQueue(subs, queue, maxSize)
    addedQueue = [subs, queue];
    if maxSize>numel(addedQueue)
        maxSize = numel(addedQueue);
    end
    if isempty(addedQueue)
       queue = [];
       return; 
    end
    [~,sortedIdx]=sort([addedQueue.mdlScore], 'descend');
    sortedQueue=addedQueue(sortedIdx);
    queue = sortedQueue(1:maxSize);
    clear sortedQueue addedQueue;
end

%> Name: removeDuplicateSubs
%>
%> Description: Given the child subs in childSubArr, this function removes
%> duplicate substructures from childSubArr, and combines the instances of
%> matching subs. If an instance can be parsed in different ways (i.e.
%> matches multiple duplicate subs), the minimum matching cost of matching is
%> saved. 
%> 
%> @param childSubArr A list of children substructures.
%>
%> @retval childSubArr The list of unique substructures.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 01.03.2015 (Converted to a standalone function, added instance augmentation)
function [childSubArr, childSubArrToEvaluate] = removeDuplicateSubs(childSubArr)
    childSubArrToEvaluate = [];
    if size(childSubArr(1).edges,1) > 1
        % Eliminate duplicate subs in final array.
        childSubEdges = cell(1, numel(childSubArr));
        for subItr = 1:numel(childSubArr)
            childSubEdges(subItr) = {sortrows(childSubArr(subItr).edges)};
        end
        subCenters = {childSubArr.centerId};
        vocabDescriptors = cellfun(@(x,y) num2str([x; y(:)]'), subCenters, childSubEdges, 'UniformOutput', false);
        [~, validChildrenIdx, IC] = unique(vocabDescriptors, 'stable');
        clear vocabDescriptors subCenters childSubEdges;
        
        % For every sub seen more than once in IC, we merge instances of
        % all matching subs.
        selectedSubs = hist(IC, numel(validChildrenIdx));
        selectedSubs = find(selectedSubs > 1);
        selectedSubBins = cell(numel(selectedSubs),1);
        for selSubItr = 1:numel(selectedSubs)
            selSubIdx = IC == selectedSubs(selSubItr);
            selectedSubBins{selSubItr} = childSubArr(selSubIdx);
        end
        updatedSubs = childSubArr(selectedSubs);
        
        parfor selSubItr = 1:numel(selectedSubs)
            % Read data for all instances.
            allInstanceCenterIdx = cat(1, selectedSubBins{selSubItr}.instanceCenterIdx);
            allInstanceChildren = cat(1, selectedSubBins{selSubItr}.instanceChildren);
            allInstanceEdges = cat(1, selectedSubBins{selSubItr}.instanceEdges);
            allInstanceMatchCosts = cat(1, selectedSubBins{selSubItr}.instanceMatchCosts);
            allInstanceCategories = cat(1, selectedSubBins{selSubItr}.instanceCategories);
            allInstanceSigns = cat(1, selectedSubBins{selSubItr}.instanceSigns);
            allInstanceDescriptors = [allInstanceCenterIdx, sort(allInstanceEdges, 2)];
            
            % Get unique instances, update minimum matching costs and write
            % everything back.
             [~, IA, IC2] = unique(allInstanceDescriptors, 'rows', 'stable');
             updatedSubs(selSubItr).instanceCenterIdx = allInstanceCenterIdx(IA);
             updatedSubs(selSubItr).instanceChildren = allInstanceChildren(IA, :);
             updatedSubs(selSubItr).instanceEdges = allInstanceEdges(IA, :);
             instanceCount = numel(IA);
             minMatchingCosts = zeros(instanceCount,1, 'single');
             % Calculate minimum matching cost for each instance.
             for instItr = 1:instanceCount
                 minMatchingCosts(instItr) = min(allInstanceMatchCosts(IC2 == instItr));
             end
            % TODO: Make matching costs optimal. Right now, they are not.
            % For the sake of efficiency, we are not looking for the minimum 
            % cost of matching for a specific instance.
%            minMatchingCosts = allInstanceMatchCosts(IA);
            updatedSubs(selSubItr).instanceMatchCosts = minMatchingCosts;
            updatedSubs(selSubItr).instanceCategories = allInstanceCategories(IA);
            updatedSubs(selSubItr).instanceSigns = allInstanceSigns(IA);
        end
        childSubArr(selectedSubs) = updatedSubs;
        
        % We return unchanged subs in childSubArr, and the subs to be
        % re-evaluated in childSubArrToEvaluate.
        if ~isempty(selectedSubs)
            childSubArrToEvaluate = childSubArr(selectedSubs);
            childSubArr = childSubArr(setdiff(IC, selectedSubs));
        end
    end
end























