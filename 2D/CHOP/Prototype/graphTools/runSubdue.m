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
function [nextVocabLevel, nextGraphLevel] = runSubdue(vocabLevel, graphLevel, nodeDistanceMatrix, edgeDistanceMatrix, categoryArrIdx, options)
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
    nsubs = options.subdue.nsubs;
    maxTime = options.subdue.maxTime;
    minSize = options.subdue.minSize;
    maxSize = options.subdue.maxSize;
    numberOfThreads = options.numberOfThreads;
    isSupervised = options.subdue.supervised;
    orgThreshold = single(options.subdue.threshold);
    regularizationParam = (options.subdue.maxSize * 2) - 1; % Maximum size of a part (n nodes + n-1 edges)
    threshold = single(options.subdue.threshold * regularizationParam) + options.singlePrecision; % Hard threshold for cost of matching two subs.
    
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
        nodeDistanceMatrix, categoryArrIdx, threshold);
 %   parentSubs = addToQueue(singleNodeSubs, parentSubs, numel(singleNodeSubs));
    parentSubs = addToQueue(singleNodeSubs, parentSubs, options.subdue.beam);
    
    %% Step 2: Main loop
    startTime = tic;
    % Let's find the adaptive threshold with the size of the minimum 
    currentSize = 2;
    minMdlScoreFinal = 0; % Preserve memory by deleting unnecessary subs.
    
    while ~isempty(parentSubs)
        adaptiveThreshold = orgThreshold * (single(currentSize) * 2 - 1) + 0.0001;
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
        
        %% All good, execution.
        display(['[SUBDUE/Parallel] Expanding subs of size ' num2str(size(parentSubs(1).edges,1)+1) '..']);
        childSubArrFinal = cell(numel(parentSubs),1);
        childSubArrToExtend = cell(numel(parentSubs),1);
        mdlScoreArrFinal = cell(numel(parentSubs),1);
        mdlScoreArrToExtend = cell(numel(parentSubs),1);
        minMdlScoreExt = 0;    % This  two will increase as more and more subs are discovered.
        setDistributions = ones(floor(numel(parentSubs)/numberOfThreads),1) * numberOfThreads;
        if rem(numel(parentSubs), numberOfThreads) > 0
            setDistributions = [setDistributions;rem(numel(parentSubs), numberOfThreads)]; %#ok<AGROW>
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
            
            % All good, continue.
            processedSet = parentSubSets{setItr};
            parfor parentItr = processedSet
                %% Step 2.2: Extend head in all possible directions into childSubs.
                childSubs = extendSub(parentSubs(parentItr), allEdges, nodeDistanceMatrix, edgeDistanceMatrix, threshold);
                if isempty(childSubs) 
                    continue;
                end
                
                % Here, we define two versions of childSubs, one for
                % extension, and one for evaluation. The evaluation version
                % has instances with a strict threshold. In order to
                % provide for the next level of discovery, we also keep 
                % one version for extension, which has instances found with
                % a more loose threshold.
                numberOfChildSubs = numel(childSubs);
                childSubsFinal = childSubs;
                validChildSubs = zeros(numberOfChildSubs,1)>0;
                for childItr = 1:numberOfChildSubs
                    validInstanceIdx = childSubsFinal(childItr).instanceMatchCosts < adaptiveThreshold;
                    if nnz(validInstanceIdx)
                        subFinal = childSubsFinal(childItr);
                        validChildSubs(childItr) = 1;
                        subFinal.instanceCenterIdx = subFinal.instanceCenterIdx(validInstanceIdx,:);
                        subFinal.instanceEdges = subFinal.instanceEdges(validInstanceIdx,:);
                        subFinal.instanceSigns = subFinal.instanceSigns(validInstanceIdx,:);
                        subFinal.instanceCategories = subFinal.instanceCategories(validInstanceIdx,:);
                        subFinal.instanceMatchCosts = subFinal.instanceMatchCosts(validInstanceIdx,:);
                        childSubsFinal(childItr) = subFinal;
                    end
                end
                childSubs = childSubs(validChildSubs);
                childSubsFinal = childSubsFinal(validChildSubs);
                
                % If no subs are remaining, continue.
                if isempty(childSubs) 
                    continue;
                end
                
                %% Step 2.3: Eliminate subs that match to more frequent subs.
                numberOfChildSubs = numel(childSubs);
                freqScores = zeros(numberOfChildSubs,1);
                for subItr = 1:numberOfChildSubs
                    freqScores(subItr) = size(childSubsFinal(subItr).instanceCenterIdx,1);
                end
                [~, sortIdx] = sort(freqScores, 'descend');
                childSubsFinal = childSubsFinal(sortIdx);
                childSubs = childSubs(sortIdx);
                validSubs = ones(numel(childSubsFinal),1)>0;
                edgeNodePairs = cat(1, childSubsFinal.edges);
                numberOfEdgesPerSub = size(childSubsFinal(1).edges,1);
                numberOfChildSubs = numel(childSubsFinal);
                edgeNodePairs = edgeNodePairs(numberOfEdgesPerSub:numberOfEdgesPerSub:(numberOfChildSubs*numberOfEdgesPerSub), :);
                for childItr = 1:(numberOfChildSubs-1)
                    if validSubs(childItr)
                        edgeNodePair1 = edgeNodePairs(childItr,:);
                        for childItr2 = (childItr+1):numberOfChildSubs
                            if validSubs(childItr2)
                                edgeNodePair2 = edgeNodePairs(childItr2,:);
                                if edgeDistanceMatrix(edgeNodePair1(1), edgeNodePair2(1)) + ...
                                        nodeDistanceMatrix(edgeNodePair1(2), edgeNodePair2(2)) < adaptiveThreshold
                                    validSubs(childItr2) = 0;
                                end
                            end
                        end
                    end
                end
                childSubs = childSubs(validSubs);
                childSubsFinal = childSubsFinal(validSubs);

                %% Step 2.4: Evaluate childSubs, find their instances.
                childSubsFinal = evaluateSubs(childSubsFinal, evalMetric, allEdges, allEdgeNodePairs, ...
                    allSigns, graphSize, overlap, mdlNodeWeight, mdlEdgeWeight, isMDLExact, isSupervised);
                
                % Assign mdl scores to childSubs, too.
                for childItr = 1:numel(childSubsFinal)
                    childSubs(childItr).mdlScore = childSubsFinal(childItr).mdlScore;
                    childSubs(childItr).normMdlScore = childSubsFinal(childItr).normMdlScore;
                end
                
                %% Sort childSubs and childSubsFinal by the mdl scores.
                mdlScores = [childSubs.mdlScore];
                [~, sortIdx] = sort(mdlScores, 'descend');
                childSubsFinal = childSubsFinal(sortIdx);
                childSubs = childSubs(sortIdx);
                mdlScores = mdlScores(sortIdx);
                
                %% Finally, eliminate subs that will not be used later on. This is done to save memory.
                mdlScoresFinal = mdlScores;
                mdlScoresExt = mdlScores;
                if currentSize == 2
                    invalidFinalSubIdx = find(mdlScoresFinal < minMdlScoreFinal, 1, 'first');
                    if ~isempty(invalidFinalSubIdx)
                        if invalidFinalSubIdx == 1
                            childSubsFinal = [];
                            mdlScoresFinal = [];
                        else
                            childSubsFinal = childSubsFinal(1, 1:(invalidFinalSubIdx-1));
                            mdlScoresFinal = mdlScoresFinal(1, 1:(invalidFinalSubIdx-1));
                        end
                    end
                    invalidExtSubIdx = find(mdlScoresExt < minMdlScoreExt, 1, 'first');
                    if ~isempty(invalidExtSubIdx)
                        if invalidExtSubIdx == 1
                            childSubs = [];
                            mdlScoresExt = [];
                        else
                            childSubs = childSubs(1, 1:(invalidExtSubIdx-1));
                            mdlScoresExt = mdlScoresExt(1, 1:(invalidExtSubIdx-1));
                        end
                    end
                end
                
                %% Save childSubs and extended subs.
                childSubArrFinal(parentItr) = {childSubsFinal};
                childSubArrToExtend(parentItr) = {childSubs};
                mdlScoreArrFinal(parentItr) = {mdlScoresFinal};
                mdlScoreArrToExtend(parentItr) = {mdlScoresExt};
            end
            
            % Update minimum limits for mdl scores for both beam and final
            % queues, so that we do not need to keep low-ranked subs in the
            % queues.
            if currentSize == 2
                allMdlScoresFinal = cat(2, mdlScoreArrFinal{:});
                if ~isempty(allMdlScoresFinal)
                    allMdlScoresFinal = sort(allMdlScoresFinal, 'descend');
                    nsubsEnd = min(numel(allMdlScoresFinal), nsubs);
                    minMdlScoreFinal = allMdlScoresFinal(nsubsEnd);
                end
                allMdlScoresExt = cat(2, mdlScoreArrToExtend{:});
                if ~isempty(allMdlScoresExt)
                    allMdlScoresExt = sort(allMdlScoresExt, 'descend');
                    beamEnd = min(numel(allMdlScoresExt), beam);
                    minMdlScoreExt = allMdlScoresExt(beamEnd);
                end
            end
        end
        
        %% Add each children group in childGroupArr into extendedSubs.
        display('[SUBDUE] Merging all children and putting them into bestSubs if they match final criteria.');
        extendedSubs = [];
        %% Step 2.4: Add childSubs to extendedSubs and bestSubs.
        childSubArrFinal = cat(2,childSubArrFinal{:});
        childSubArrToExtend = cat(2, childSubArrToExtend{:});
        % Add children to the both queues.
        if ~isempty(childSubArrFinal)
            % Remove duplicate nodes from the children subs.
            if size(childSubArrFinal(1).edges,1) > 1
                % Eliminate duplicate subs in final array.
                childSubEdges = cell(1, numel(childSubArrFinal));
                for subItr = 1:numel(childSubArrFinal)
                    childSubEdges(subItr) = {sortrows(childSubArrFinal(subItr).edges)};
                end
                subCenters = {childSubArrFinal.centerId};
                vocabDescriptors = cellfun(@(x,y) num2str([x; y(:)]'), subCenters, childSubEdges, 'UniformOutput', false);
                [~, validChildrenIdx, ~] = unique(vocabDescriptors, 'stable');
                clear vocabDescriptors subCenters childSubEdges;
                childSubArrFinal = childSubArrFinal(validChildrenIdx);
                childSubArrToExtend = childSubArrToExtend(validChildrenIdx);
            end
        
            % Check for maxSize to put to extendedSubs (parentSubs for next level).
            if currentSize < maxSize
                extendedSubs = addToQueue(childSubArrToExtend, extendedSubs, beam);
            end
            % Check for minSize to put to final sub array.
            if currentSize >= minSize
                bestSubs = addToQueue(childSubArrFinal, bestSubs, nsubs); 
            end
        end
        clear childSubArrFinal childSubArrToExtend;
        
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
%        remainingBestSubs = 1:numberOfBestSubs;
%        IA = 1:numberOfInstances;
        %% First, we eliminate the instances having exactly the same node list.
        % Learn the maximum children count in instances.
        maxSubSize = 1;
        for bestSubItr = 1:numberOfBestSubs
            maxSubSize = max(maxSubSize, (size(bestSubs(bestSubItr).edges,1)+1));
        end
        
        % Find the children for each instance, and use it as a descriptor.
        % We'll then find unique descriptors, which is the initial phase of
        % inhibition.
        remainingInstanceLabels = zeros(numberOfInstances, 1);
        instanceChildrenDescriptors = zeros(numberOfInstances,maxSubSize, 'int32');
       instanceOffset = 1;
       for bestSubItr = 1:numberOfBestSubs
           centerIdx = bestSubs(bestSubItr).instanceCenterIdx;
           centerIdxCellArr = num2cell(centerIdx);
           edges = {allEdges(centerIdx).adjInfo}';
           edgeIdx = bestSubs(bestSubItr).instanceEdges;
           edgeIdx = mat2cell(edgeIdx, ones(size(edgeIdx,1),1), size(edgeIdx,2));
           instanceEdges = cellfun(@(x,y) x(y,:), edges, edgeIdx, 'UniformOutput', false);
           instanceChildren = cellfun(@(x,y) sort([x, y(:,2)']), centerIdxCellArr, instanceEdges, 'UniformOutput', false);
           instanceChildren = cat(1, instanceChildren{:});
           instanceEndOffset = instanceOffset + numel(centerIdx) - 1;
           instanceChildrenDescriptors(instanceOffset:instanceEndOffset, 1:size(instanceChildren,2)) = instanceChildren;
           remainingInstanceLabels(instanceOffset:instanceEndOffset) = bestSubItr;
           instanceOffset = instanceEndOffset + 1;
       end
       % Find unique rows, which correspond to unique instances.
       [~, IA, ~] = unique(instanceChildrenDescriptors, 'rows', 'stable');
       
       % Get the remaining sub count.
       remainingBestSubs = unique(remainingInstanceLabels(IA))';
       numberOfBestSubs = numel(remainingBestSubs);
       numberOfInstances = numel(IA);
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
               centerIdxCellArr = num2cell(centerIdx);
               numberOfInstances = numel(centerIdx);
               edges = {allEdges(centerIdx).adjInfo}';
               edgeIdx = bestSubs(bestSubItr).instanceEdges(remainingInstanceIdx,:);
               edgeIdx = mat2cell(edgeIdx, ones(numberOfInstances, 1), size(edgeIdx,2));
               instanceEdges = cellfun(@(x,y) x(y,:), edges, edgeIdx, 'UniformOutput', false);
     %          clear edges edgeIdx;
               instanceChildren = cellfun(@(x,y) [x, y(:,2)'], centerIdxCellArr, instanceEdges, 'UniformOutput', false);
      %         clear centerIdxCellArr;
               instanceSigns = num2cell(allSigns(centerIdx));

               % Assign fields to graphs.
               [graphLevel(actualInstanceOffset:(actualInstanceOffset + numberOfInstances-1)).labelId] = deal(labelIds{:});
               [graphLevel(actualInstanceOffset:(actualInstanceOffset + numberOfInstances-1)).children] = deal(instanceChildren{:});
               [graphLevel(actualInstanceOffset:(actualInstanceOffset + numberOfInstances-1)).sign] = deal(instanceSigns{:});
     %          clear labelIds instanceChildren instanceSigns centerIdx;
               instanceOffset = instanceEndOffset + 1;
               actualInstanceOffset = actualInstanceOffset + numberOfInstances;
               remainingSubItr = remainingSubItr + 1;
           end
        end
    else
        vocabLevel = [];
        graphLevel = [];
    end
   
    
   clear allEdges allEdgeNodePairs allSigns bestSubs;
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
function [extendedSubs] = extendSub(sub, allEdges, nodeDistanceMatrix, edgeDistanceMatrix, threshold)
    centerIdx = sub.instanceCenterIdx;
    subAllEdges = {allEdges(centerIdx).adjInfo}';
    % Get unused edges from sub's instances.
    allUsedEdgeIdx = sub.instanceEdges;
    if ~isempty(allUsedEdgeIdx)
        allUsedEdgeIdx = mat2cell(allUsedEdgeIdx, ones(size(allUsedEdgeIdx,1),1), size(allUsedEdgeIdx,2));
        allUnusedEdges = cellfun(@(x,y) x(setdiff(1:size(x,1),y),:), subAllEdges, allUsedEdgeIdx, 'UniformOutput', false);
    else
        allUnusedEdges = subAllEdges;
    end
 %   clear allUsedEdgeIdx subAllEdges;
    
    % Record which edge belongs to which instance. 
    allEdgeInstanceIds = zeros(sum(cellfun(@(x) size(x,1), allUnusedEdges)),1);
    allEdgePrevCosts = zeros(size(allEdgeInstanceIds), 'single');
    itrOffset = 1;
    unusedEdgeCount = cellfun(@(x) size(x,1), allUnusedEdges);
    for itr = 1:numel(allUnusedEdges)
        beginOffset = itrOffset;
        endOffset = (beginOffset+(unusedEdgeCount(itr)-1));
        allEdgeInstanceIds(beginOffset:endOffset) = itr;
        allEdgePrevCosts(beginOffset:endOffset) = sub.instanceMatchCosts(itr);
        itrOffset = itrOffset + unusedEdgeCount(itr);
    end         
    allUnusedEdges = cat(1, allUnusedEdges{:});
    
    % If no edges remain, exit.
    if isempty(allUnusedEdges)
        extendedSubs = [];
        return; 
    end
    
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
        edgesToExtend = allUnusedEdges(edgesToExtendIdx,:);
        edgeInstanceIds = allEdgeInstanceIds(edgesToExtendIdx);
        
        % Assign fields related to the new instances.
        newSub.instanceCenterIdx = sub.instanceCenterIdx(edgeInstanceIds);
        newSub.instanceSigns = sub.instanceSigns(edgeInstanceIds);
        newSub.instanceCategories = sub.instanceCategories(edgeInstanceIds);
        newSub.instanceMatchCosts = edgesToExtendCosts(edgesToExtendIdx);
        existingEdges = sub.instanceEdges;
        combinedEdges = zeros(numel(edgeInstanceIds), size(existingEdges,2)+1, 'uint8');
        if ~isempty(existingEdges)
            existingEdges = existingEdges(edgeInstanceIds);
            combinedEdges(:, 1:size(existingEdges,2)) =  existingEdges;
        end
        
        % Each added edge actually means a new instance. 
        for instanceItr = 1:numel(edgeInstanceIds)
            instanceEdges = allEdges(centerIdx(edgeInstanceIds(instanceItr))).adjInfo;
            if isempty(instanceEdges)
               continue; 
            end
            addedEdgeIdx = find(instanceEdges(:, 2) == edgesToExtend(instanceItr,2));
            combinedEdges(instanceItr, end) = addedEdgeIdx;
        end
        newSub.instanceEdges = combinedEdges;
        
        % All instances assigned, good to go.
        extendedSubs(edgeTypeItr) = newSub;
 %       clear newInstances newSub edgesToExtendCosts edgesToExtendIdx edgesToExtend edgeInstanceIds instanceEdges;
    end
 %   clear subAllEdges;
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
    for subItr = 1:numberOfSubs
        % Find the weight of this node, by taking the max of the category distribution. 
        if isSupervised
            categoryArr = double(subs(subItr).instanceCategories);
            weight = nnz(categoryArr == mode(categoryArr)) / numel(categoryArr);
        else
            weight = 1;
        end
        
        % We compress the object graph using the children, and the
        % edges they are involved. 
        subScore = weight * getSubScore(subs(subItr), allEdges, allEdgeNodePairs, evalMetric, ...
           allSigns, mdlNodeWeight, mdlEdgeWeight, ....
            overlap, isMDLExact);

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